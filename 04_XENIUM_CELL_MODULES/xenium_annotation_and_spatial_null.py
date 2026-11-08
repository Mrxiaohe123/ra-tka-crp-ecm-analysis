# PUBLIC RELEASE CANDIDATE COPY
# Source: STAGE7Q_SUPPLEMENTARY_EVIDENCE_PATCH/SOURCE_CODE/08_stage7d_qa_spatial.py
# Source SHA256: e35ebb1b3ea32a5a7c85b88d403463f42730531d7610a03b6dac6870f028bfa3
# Purpose: Xenium annotation sensitivity and fixed spatial null model
# This copy does not overwrite the source and contains no study input data.

#!/usr/bin/env python3
"""Stage 7D single_donor_pilot annotation sensitivity, ECM source and spatial null model."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

ALT_MARKERS = {
    "Fibroblast": ["PDPN", "THY1", "PRG4", "CXCL12", "FAP", "CD34"],
    "Myeloid": ["C1QA", "C1QC", "CD68", "FCGR3A", "TYROBP", "AIF1", "MS4A7"],
    "T_cell": ["TRAC", "IL7R", "LTB", "CD8A"],
    "B_cell": ["MS4A1", "CD79A"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("h5ad", type=Path); ap.add_argument("cells", type=Path)
    ap.add_argument("annotation", type=Path); ap.add_argument("coverage", type=Path)
    ap.add_argument("out", type=Path); ap.add_argument("--seed", type=int, default=20261029)
    args=ap.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    a=ad.read_h5ad(args.h5ad); cells=pd.read_csv(args.cells); old=pd.read_csv(args.annotation)
    cells["cell_id"]=cells["cell_id"].astype(str); old["cell_id"]=old["cell_id"].astype(str)
    obs_names=a.obs_names.astype(str).tolist()
    if set(obs_names)!=set(cells.cell_id) or set(obs_names)!=set(old.cell_id): raise SystemExit('ID sets do not match')
    cells=cells.set_index('cell_id').loc[obs_names]; old=old.set_index('cell_id').loc[obs_names]
    X=a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
    X=X.astype(float,copy=False); totals=X.sum(axis=1); detected=(X>0).sum(axis=1)
    norm=np.log1p(X/np.maximum(totals[:,None],1)*1e4); var=np.asarray(a.var_names.astype(str)); idx={g:i for i,g in enumerate(var)}
    cov=pd.read_csv(args.coverage); frozen=cov.gene.astype(str).tolist(); frozen_set=set(frozen)

    score={}; present={}
    for lab,genes in ALT_MARKERS.items():
        genes=[g for g in genes if g not in frozen_set]
        p=[g for g in genes if g in idx]; present[lab]=p
        score[lab]=norm[:,[idx[g] for g in p]].mean(axis=1) if p else np.zeros(len(norm))
    s=pd.DataFrame(score,index=old.index); arr=np.sort(s.to_numpy(),axis=1); top=s.idxmax(axis=1); margin=arr[:,-1]-arr[:,-2]
    nd=np.zeros(len(norm),dtype=int)
    for lab,p in present.items(): nd+=(X[:,[idx[g] for g in p]]>0).sum(axis=1) if p else 0
    alt=np.where((nd>=2)&(margin>=0.05),top,"Unassigned")
    comp=pd.DataFrame({
        'cell_id':old.index,'original_label':old.annotation.to_numpy(),'alternative_label':alt,
        'agreement':old.annotation.to_numpy()==alt,'original_confidence':old.annotation_confidence.to_numpy(),
        'total_counts':totals,'n_genes_detected':detected,'alternative_marker_detected_n':nd,
        'alternative_top_score':arr[:,-1],'alternative_score_margin':margin,
    })
    comp.to_csv(args.out/'SINGLE_DONOR_PILOT_ANNOTATION_COMPARISON.csv',index=False)

    # ECM source table by exact Stage 2 frozen category and both label schemes.
    rows=[]
    for scheme, labels in [('original_marker',old.annotation.to_numpy()),('alternative_non_ecm_marker',alt)]:
        for cat in cov.frozen_category.unique():
            genes=cov.loc[cov.frozen_category==cat,'gene'].astype(str).tolist(); p=[g for g in genes if g in idx]
            if not p: continue
            vals=norm[:,[idx[g] for g in p]]; gscore=vals.mean(axis=1)
            for lab in sorted(set(labels)):
                take=labels==lab; vv=vals[take]
                rows.append({'annotation_scheme':scheme,'label':lab,'module':cat,'n_cells':int(take.sum()),'n_genes_on_panel':len(p),'mean_module_score':float(gscore[take].mean()),'median_module_score':float(np.median(gscore[take])),'mean_gene_detection_rate':float((vv>0).mean()),'median_genes_detected_per_cell':float(np.median((vv>0).sum(axis=1)))})
    pd.DataFrame(rows).to_csv(args.out/'ECM_CELL_SOURCE_REASSESSMENT.csv',index=False)

    # Spatial null: fixed 4x4 equal-width bins, permute labels inside each bin.
    labels=old.annotation.to_numpy().astype(object); xy=cells[['x_centroid','y_centroid']].to_numpy(float)
    fib=labels=='Fibroblast'; mye=labels=='Myeloid'; radius=100.0
    def metrics(lab):
        fi=np.where(lab=='Fibroblast')[0]; mi=np.where(lab=='Myeloid')[0]
        if len(fi)==0 or len(mi)==0:return np.nan,np.nan
        d,_=cKDTree(xy[mi]).query(xy[fi],k=1)
        return float(np.median(d)),float(np.mean(d<=radius))
    obs_d,obs_f=metrics(labels)
    xmin,ymin=xy.min(axis=0); xmax,ymax=xy.max(axis=0)
    bx=np.minimum(3,np.floor((xy[:,0]-xmin)/max(xmax-xmin,1e-9)*4).astype(int)); by=np.minimum(3,np.floor((xy[:,1]-ymin)/max(ymax-ymin,1e-9)*4).astype(int)); bins=by*4+bx
    rng=np.random.default_rng(args.seed); ndraw=1000; null=[]
    for _ in range(ndraw):
        perm=labels.copy()
        for b in range(16):
            ii=np.where(bins==b)[0]
            if len(ii)>1: perm[ii]=perm[rng.permutation(ii)]
        null.append(metrics(perm))
    null=np.asarray(null); out=[]
    for metric,ov,vals,direction in [('nearest_distance_median_um',obs_d,null[:,0],'low'),('fraction_fibro_with_myeloid_within_100um',obs_f,null[:,1],'high')]:
        out.append({'metric':metric,'observed':ov,'null_mean':float(np.nanmean(vals)),'null_q025':float(np.nanquantile(vals,.025)),'null_q975':float(np.nanquantile(vals,.975)),'n_permutations':ndraw,'seed':args.seed,'p_permutation_one_sided':float((1+np.sum(vals<=ov if direction=='low' else vals>=ov))/(ndraw+1)),'null_model':'fixed 4x4 spatial-bin label permutation preserving within-bin label counts'})
    pd.DataFrame(out).to_csv(args.out/'ECM_SPATIAL_NULL_MODEL_RESULTS.csv',index=False)
    (args.out/'stage7d_summary.json').write_text(json.dumps({'n_cells':len(labels),'alternative_marker_panel_present':present,'annotation_counts_original':pd.Series(labels).value_counts().to_dict(),'annotation_counts_alternative':pd.Series(alt).value_counts().to_dict(),'n_agreement':int(comp.agreement.sum()),'agreement_fraction':float(comp.agreement.mean()),'unassigned_original_n':int((labels=='Unassigned').sum()),'unassigned_original_median_counts':float(np.median(totals[labels=='Unassigned'])),'unassigned_original_median_genes':float(np.median(detected[labels=='Unassigned'])),'spatial_observed':{'nearest_distance_median_um':obs_d,'fraction_within_100um':obs_f},'spatial_null':out},indent=2),encoding='utf-8')

if __name__=='__main__': main()
