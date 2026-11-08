# PUBLIC RELEASE CANDIDATE COPY
# Source: VERSION_C_FINAL_CLINICAL_SOURCE_ADJUDICATION/reproduce_source_and_auxiliary_MAR.py
# Source SHA256: 617ee3f84dd93c0451c1b45674ce217c27afa8bd4a839950e0463524e3fb793a
# Purpose: Auxiliary MAR/MNAR source audit; institution-controlled workbook required
# This copy does not overwrite the source and contains no study input data.

#!/usr/bin/env python3
"""Version C: final-workbook aggregate audit and auxiliary MAR multiple imputation.

Requires artifact_tool, numpy, scipy and statsmodels. Reads a local Excel workbook;
exports only aggregated results, never individual patient records or ID joins.
Run: python reproduce_source_and_auxiliary_MAR.py --workbook PATH --output-dir DIR
Optional: --prior-snapshot ZIP (historical V2 supplement source)

The auxiliary MAR method is a Bayesian multivariate-normal imputation of log CRP
(M3, M6, M12) conditional on log pre-CRP, age, sex, RA duration, BMI and
baseline IL-6 inhibitor use. It incorporates each available intermediate CRP.
4 independent chains × 250 imputations; Rubin pooling of paired log M12/Pre.
One documented death before M12 is excluded (conditional N=79 analysis).
All other missing outcomes are imputed conditional on assumptions, WITHOUT
asserting that the 29 persons with unknown missingness reasons survived.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, math, os, warnings, zipfile
from collections import Counter
from pathlib import Path
from datetime import date, timedelta
import numpy as np
from scipy.stats import invwishart, invgamma, t as student_t
from artifact_tool import Blob, SpreadsheetFile
import statsmodels.api as sm

CRP_VISITS = ['Pre', 'M3', 'M6', 'M12']
CONTROLS = ['Age','Sex','RA_Duration_y','BMI','Pre_CRP','Base_IL6i']
SEED_BASE = 20261106

def read_records(path):
    book=SpreadsheetFile.import_xlsx(Blob.load(str(path)))
    source=book.worksheets.get_item_at(0)
    raw=source.get_range('A2:AW82').values
    fields=raw[0]
    records=[dict(zip(fields,line)) for line in raw[1:] if str(line[0] or '').startswith('RA_TKA_')]
    assert len(records)==80 and len({r['Patient_ID'] for r in records})==80
    return records,fields

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def cvalue(v):
    if v is None or v=='':return None
    return round(float(v),8)

def write_csv(path,head,rows):
    with open(path,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(head);w.writerows(rows)

def clinical_summaries(rs):
    observed={t:sum(x[t+'_CRP'] is not None for x in rs) for t in CRP_VISITS}
    flags={'sex':Counter(x['Sex'] for x in rs),'laterality':Counter(x['Bilateral'] for x in rs),
           'il6i':Counter(x['Base_IL6i'] for x in rs),'exclude':Counter(x['Exclude'] for x in rs)}
    assert observed=={'Pre':80,'M3':75,'M6':57,'M12':46},observed
    death=[r for r in rs if r['Exclude']==1 and ('死亡' in str(r['Exclude_Reason']))]
    assert len(death)==1 and death[0]['M12_CRP'] is None
    bilats=[r for r in rs if r['Bilateral']==2]
    gaps=[r['Surg_Date']-r['Pre_Date'] for r in bilats]
    return observed, flags, death[0], gaps

def mixed_model(rs):
    yy=[]; xx=[];groups=[]
    for i,r in enumerate(rs):
        for j,t in enumerate(CRP_VISITS):
            crp=r[t+'_CRP']
            if isinstance(crp,(int,float)) and crp>0:
                yy.append(np.log(float(crp)))
                xx.append([1.0]+[float(j==k) for k in (1,2,3)])
                groups.append(i)
    yy=np.array(yy); xx=np.array(xx)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res=sm.MixedLM(yy,xx,groups=np.array(groups)).fit(reml=True,method='bfgs',disp=False,maxiter=500)
    coeff=np.asarray(res.fe_params)
    se=np.asarray(res.bse_fe)
    z=1.959963984540054
    ci=np.exp(np.column_stack([coeff[1:]-z*se[1:], coeff[1:]+z*se[1:]]))
    return len(rs),len(yy),np.exp(coeff[1:]),ci,res.converged

def mixed_model_adjusted(rs):
    yy=[];xx=[];groups=[]
    for i,r in enumerate(rs):
        for j,t in enumerate(CRP_VISITS):
            crp=r[t+'_CRP']
            if isinstance(crp,(int,float)) and crp>0:
                yy.append(np.log(float(crp)))
                xx.append([1.0]+[float(j==q) for q in (1,2,3)]+[
                    float(r['Age']),float(r['Sex']==2),float(r['RA_Duration_y']),float(r['BMI']),float(r['Base_IL6i'])])
                groups.append(i)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        result=sm.MixedLM(np.array(yy),np.array(xx),groups=np.array(groups)).fit(reml=True,method='bfgs',disp=False,maxiter=500)
    b=np.asarray(result.fe_params)[1:4]
    se=np.asarray(result.bse_fe)[1:4]
    z=1.959963984540054
    return len(rs),len(yy),np.exp(b),np.exp(np.column_stack([b-z*se,b+z*se])),result.converged

def compare_snapshot(rs,prior_zip):
    final=Counter(tuple(cvalue(r[t+'_CRP']) for t in CRP_VISITS) for r in rs)
    notes=[]
    if prior_zip and Path(prior_zip).exists():
        with zipfile.ZipFile(prior_zip) as f:
            candidates=[s for s in f.namelist() if s.endswith('s1_observation_matrix_80x4.csv')]
            if candidates:
                with f.open(candidates[0]) as g:
                    old=list(csv.DictReader(io.TextIOWrapper(g,encoding='utf-8-sig')))
                earlier=Counter(tuple(cvalue(r[t+'_recorded_value']) for t in CRP_VISITS) for r in old)
                notes.append(('Historical V2 supplement snapshot',len(old),sum((final&earlier).values()),sum((final-earlier).values()),sum((earlier-final).values()),'Matched full 4-visit CRP tuples as multisets; anonymized plot IDs are not a patient ID crosswalk'))
                return notes
    notes.append(('Historical complete workbook',None,None,None,None,'NOT_AVAILABLE; full cell-by-cell comparison cannot be performed'))
    return notes

def construct_x(rs):
    n=len(rs)
    logpre=np.array([np.log(float(r['Pre_CRP'])) for r in rs])
    # Sex: 1=male, 2=female; standardize centered predictors using sample SD.
    raw=np.column_stack([
        logpre,
        [float(r['Age']) for r in rs],
        [float(r['Sex']==2) for r in rs],
        [float(r['RA_Duration_y']) for r in rs],
        [float(r['BMI']) for r in rs],
        [float(r['Base_IL6i']) for r in rs],
    ])
    means=np.mean(raw,axis=0);std=np.std(raw,axis=0,ddof=1)
    assert np.all(std>0)
    x=np.column_stack([np.ones(n),(raw-means)/std])
    y=np.array([[np.log(float(r[t+'_CRP'])) if r[t+'_CRP'] is not None else np.nan for t in ['M3','M6','M12']] for r in rs])
    assert np.all(np.isfinite(y[~np.isnan(y)]))
    return x,y,logpre,means,std

def single_chain(x,obs_y, logpre,seed, draws=250,burn=550,thin=3):
    rng=np.random.default_rng(seed)
    n,k=obs_y.shape;p=x.shape[1]
    seen=np.isfinite(obs_y)
    y=obs_y.copy()
    col_means=np.nanmean(obs_y,axis=0)
    for j in range(k):y[~seen[:,j],j]=col_means[j]
    # Weakly regularizing normal-inverse-Wishart prior on Bayesian MVN imputation.
    b0=np.zeros((p,k)); b0[0,:]=col_means
    v0inv=np.diag([0.005]+[0.05]*(p-1))
    v=np.linalg.inv(x.T@x+v0inv)
    lv=np.linalg.cholesky(v)
    df0=k+3
    psi0=np.eye(k)*0.6
    groups=[]
    for pattern in set(tuple(row) for row in seen.tolist()):
        inds=np.where(np.all(seen==np.array(pattern)[None,:],axis=1))[0]
        mis=np.where(~np.array(pattern))[0]
        if len(mis):groups.append((inds,np.where(np.array(pattern))[0],mis))
    out=[]
    total=burn+draws*thin
    for step in range(total):
        mn=v@(x.T@y+v0inv@b0)
        resid=y-x@mn
        db=mn-b0
        pscale=psi0+resid.T@resid+db.T@v0inv@db
        sigma=invwishart.rvs(df=df0+n,scale=pscale,random_state=rng)
        chol=np.linalg.cholesky(sigma)
        b=mn+lv@rng.standard_normal((p,k))@chol.T
        mu=x@b
        for inds, ob, mi in groups:
            if len(ob):
                soo=sigma[np.ix_(ob,ob)]
                smo=sigma[np.ix_(mi,ob)]
                a=smo@np.linalg.inv(soo)
                center=mu[np.ix_(inds,mi)]+(y[np.ix_(inds,ob)]-mu[np.ix_(inds,ob)])@a.T
                cov=sigma[np.ix_(mi,mi)]-a@sigma[np.ix_(ob,mi)]
            else:
                center=mu[np.ix_(inds,mi)]
                cov=sigma[np.ix_(mi,mi)]
            cov=(cov+cov.T)/2
            value=rng.multivariate_normal(np.zeros(len(mi)),cov,size=len(inds))
            y[np.ix_(inds,mi)]=center+value
        if step>=burn and (step-burn)%thin==0:
            delta=y[:,2]-logpre
            qhat=delta.mean()
            uhat=delta.var(ddof=1)/n
            out.append((qhat,uhat))
    return np.array(out)

def baseline_only_draws(x,obs_y,logpre,seed,draws=250):
    """Conjugate Bayesian linear predictive MI for missing log M12 CRP only.
    Baseline predictors identical to full auxiliary model; no M3/M6 inputs.
    """
    rng=np.random.default_rng(seed)
    observed=np.isfinite(obs_y)
    xo=x[observed,:]
    yo=obs_y[observed]
    xm=x[~observed,:]
    n,p=x.shape
    beta0=np.zeros(p); beta0[0]=np.mean(yo)
    v0inv=np.diag([0.005]+[0.05]*(p-1))
    v=np.linalg.inv(xo.T@xo+v0inv)
    mn=v@(xo.T@yo+v0inv@beta0)
    alpha=2.0+len(yo)/2
    scale=0.3+0.5*(np.sum((yo-xo@mn)**2)+(mn-beta0)@v0inv@(mn-beta0))
    lv=np.linalg.cholesky(v)
    out=[]
    for _ in range(draws):
        s2=invgamma.rvs(a=alpha,scale=scale,random_state=rng)
        beta=mn+np.sqrt(s2)*lv@rng.standard_normal(p)
        completed=obs_y.copy()
        completed[~observed]=xm@beta+np.sqrt(s2)*rng.standard_normal(len(xm))
        d=completed-logpre
        out.append((d.mean(),d.var(ddof=1)/n))
    return np.array(out)

def pool_rubin(est,within,dfcomplete):
    m=len(est); qbar=float(np.mean(est));ubar=float(np.mean(within));b=float(np.var(est,ddof=1));tt=ubar+(1+1/m)*b
    lam=(1+1/m)*b/tt
    if lam>0:
        old=(m-1)/lam**2
        obs=((dfcomplete+1)/(dfcomplete+3))*dfcomplete*(1-lam)
        dfr=1/(1/old+1/obs)
    else:dfr=dfcomplete
    qcrit=float(student_t.ppf(.975,dfr));lo=qbar-qcrit*np.sqrt(tt);hi=qbar+qcrit*np.sqrt(tt)
    return {'n_imputations':m,'log_gmr':qbar,'gmr':float(np.exp(qbar)),'ci_low':float(np.exp(lo)),'ci_high':float(np.exp(hi)),
    'se_log':float(np.sqrt(tt)),'df':float(dfr),'within_variance':ubar,'between_variance':b,'fraction_missing_information':float(lam),
    'p_two_sided':float(2*student_t.sf(abs(qbar/np.sqrt(tt)),dfr))}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--workbook', type=Path, required=True)
    ap.add_argument('--prior-snapshot', type=Path, required=False, default=None)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--draws-per-chain',type=int,default=250)
    ap.add_argument('--burn',type=int,default=550)
    ap.add_argument('--thin',type=int,default=3)
    args=ap.parse_args()
    dest=Path(args.output_dir);dest.mkdir(parents=True,exist_ok=True)
    rec,_=read_records(args.workbook)
    observed, flags,death,gaps=clinical_summaries(rec)
    write_csv(dest/'FINAL_WORKBOOK_FIELD_AUDIT.csv',['field','available_n','missing_n','source_status'],
              [(k,sum(r.get(k) is not None for r in rec),sum(r.get(k) is None for r in rec),'FINAL_WORKBOOK_ONLY; earlier full workbook unavailable')
               for k in CONTROLS+['Surg_Date','Pre_Date','Bilateral','Pre_VAS','M6_VAS','M3_CRP','M6_CRP','M12_CRP','Exclude','Exclude_Reason']])
    write_csv(dest/'HISTORICAL_SNAPSHOT_CRP_MATCH.csv',['source','historical_n','full_tuple_match_n','new_unmatched_n','historical_unmatched_n','interpretation'],compare_snapshot(rec,args.prior_snapshot))
    write_csv(dest/'LATERALITY_SOURCE_DECISION.csv',['code','meaning_confirmed_by_author','patients_n','percent','validation'],
              [(0,'Right unilateral TKA',flags['laterality'][0],round(100*flags['laterality'][0]/80,2),'Final workbook aggregate'),
               (1,'Left unilateral TKA',flags['laterality'][1],round(100*flags['laterality'][1]/80,2),'Final workbook aggregate'),
               (2,'Bilateral TKA',flags['laterality'][2],round(100*flags['laterality'][2]/80,2),'Final workbook aggregate')])
    res=[]
    for label,rows in [('N80_MAIN',rec),('N72_EXCLUDE_IL6I',[r for r in rec if r['Base_IL6i']==0]),('N75_EXCLUDE_FLAG5',[r for r in rec if r['Exclude']==0])]:
        n,nobs,gmr,ci,converged=mixed_model(rows)
        for j,t in enumerate(['M3','M6','M12']):
            res.append((label,n,nobs,t,round(float(gmr[j]),7),round(float(ci[j,0]),7),round(float(ci[j,1]),7),converged))
    n,nobs,gmr,ci,converged=mixed_model_adjusted(rec)
    for j,t in enumerate(['M3','M6','M12']):
        res.append(('N80_BASELINE_ADJUSTED',n,nobs,t,round(float(gmr[j]),7),round(float(ci[j,0]),7),round(float(ci[j,1]),7),converged))
    write_csv(dest/'MIXED_MODEL_FINAL_WORKBOOK_REPRODUCTION.csv',['analysis','patients_n','total_observations','contrast','gmr','95ci_lower','95ci_upper','model_converged'],res)
    death_index=rec.index(death)
    survivors_cond=[r for r in rec if r is not death]
    assert len(survivors_cond)==79
    x,y,lp,means,std=construct_x(survivors_cond)
    # actual observed counts and auxiliary distributions, stratified by M12 availability
    aux=[]
    for has12 in (True,False):
        rr=[r for r in rec if (r['M12_CRP'] is not None)==has12]
        for time in ['M3','M6']:
            a=np.array([r[time+'_CRP'] for r in rr if r[time+'_CRP'] is not None],dtype=float)
            aux.append(('M12 observed' if has12 else 'M12 missing',len(rr),time,len(a),
                        round(float(np.median(a)),4),round(float(np.quantile(a,.25)),4),round(float(np.quantile(a,.75)),4)))
    write_csv(dest/'M12_INTERMEDIATE_OBSERVATION_DIAGNOSTICS.csv',['M12_status','group_n','visit','observed_n','crp_median_mgL','crp_q1_mgL','crp_q3_mgL'],aux)
    # Proper Bayesian multivariate joint-model imputation. 4 independent chains.
    chains=[]
    for chain in range(4):
        vals=single_chain(x,y,lp,seed=SEED_BASE+chain*101,draws=args.draws_per_chain,burn=args.burn,thin=args.thin)
        chains.append(vals)
    vals=np.concatenate(chains,axis=0)
    result=pool_rubin(vals[:,0],vals[:,1],dfcomplete=len(survivors_cond)-1)
    write_csv(dest/'M12_AUXILIARY_MAR_RESULTS.csv',
              ['analysis','conditional_n','M12_observed_n','M12_imputed_n','M3_missing_n','M6_missing_n','imputations','gmr','ci_lower','ci_upper','p_two_sided','log_estimate','pooled_se_log','df','within_var','between_var','fmi'],
              [('JOINT_LOGCRP_M3_M6_M12_MAR',len(survivors_cond),int(np.isfinite(y[:,2]).sum()),int(np.isnan(y[:,2]).sum()),int(np.isnan(y[:,0]).sum()),int(np.isnan(y[:,1]).sum()),len(vals),
                result['gmr'],result['ci_low'],result['ci_high'],result['p_two_sided'],result['log_gmr'],result['se_log'],result['df'],result['within_variance'],result['between_variance'],result['fraction_missing_information'])])
    base_chains=[baseline_only_draws(x,y[:,2],lp,seed=SEED_BASE+5000+ci*101,draws=args.draws_per_chain) for ci in range(4)]
    basevals=np.concatenate(base_chains,axis=0)
    base_result=pool_rubin(basevals[:,0],basevals[:,1],dfcomplete=len(survivors_cond)-1)
    write_csv(dest/'M12_MAR_BASELINE_VS_INTERMEDIATE_COMPARISON.csv',
              ['method','conditional_n','M12_imputed_n','imputations','gmr','95ci_lower','95ci_upper','p_two_sided','log_gmr','se_log','assumption'],
              [('BASELINE_ONLY_CONJUGATE_MI',len(survivors_cond),int(np.isnan(y[:,2]).sum()),len(basevals),base_result['gmr'],base_result['ci_low'],base_result['ci_high'],base_result['p_two_sided'],base_result['log_gmr'],base_result['se_log'],
                'M12 MAR conditional on baseline clinical covariates; no intermediate CRP'),
               ('AUXILIARY_JOINT_M3_M6_M12_MI',len(survivors_cond),int(np.isnan(y[:,2]).sum()),len(vals),result['gmr'],result['ci_low'],result['ci_high'],result['p_two_sided'],result['log_gmr'],result['se_log'],
                'M12 MAR conditional on baseline covariates and observed M3/M6 CRP, joint multivariate normal')])
    chain_stats=[(i+1,float(np.exp(np.mean(c[:,0]))),float(np.mean(c[:,0])),float(np.std(c[:,0],ddof=1))) for i,c in enumerate(chains)]
    write_csv(dest/'M12_AUXILIARY_MAR_CHAIN_DIAGNOSTICS.csv',['chain','mean_gmr','mean_log_gmr','sd_imputed_log_estimate'],chain_stats)
    write_csv(dest/'M12_AUXILIARY_MAR_REPRO_SETTINGS.csv',['setting','value'],[
        ('imputation_model','Bayesian multivariate normal for log M3, log M6, log M12 CRP conditional on baseline covariates'),
        ('baseline_covariates','log(Pre CRP), age, sex, RA duration, BMI, baseline IL6i'),
        ('M3_and_M6_role','observed intermediate CRP values enter joint imputation; missing intermediate values are imputed'),
        ('M12_target','geometric mean paired M12/Pre across N79 conditional records, excluding documented death'),
        ('priors','weak normal-inverse-Wishart: V0inv=diag(.005,.05×6), nu0=6, Psi0=0.6 I'),
        ('chain_count',4),('draws_per_chain',args.draws_per_chain),('burn',args.burn),('thin',args.thin),('seed_base',SEED_BASE),
        ('pooling','Rubin total variance; Barnard-Rubin df; complete-data pairwise mean-log-difference within SE'),
        ('M12_unobserved_survival','29 patients missing cause not verified alive; model conditional/assumption-based, not survivor-cohort estimand')])
    qa={'filename':str(args.workbook),'sha256':sha(args.workbook),'n':80,'crp_available':observed,
        'n_bilateral':5,'bilateral_preoperative_date_gap_days':sorted(gaps),
        'n_il6i':flags['il6i'][1],'n_exclude5':flags['exclude'][1],
        'death_documented':1,'death_after_m6_no_m12':True,
        'n_imputation_target':79,'n_imputed_m12':int(np.isnan(y[:,2]).sum()),
        'imputation_result':result,'baseline_only_MAR_comparison':base_result,'per_chain_mean_gmr':chain_stats,
        'full_historical_workbook_available':False,'historical_crp_snapshot_available':True}
    with open(dest/'QA_AND_METHODS.json','w',encoding='utf-8') as f:json.dump(qa,f,ensure_ascii=False,indent=2)
    print(json.dumps(qa,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
