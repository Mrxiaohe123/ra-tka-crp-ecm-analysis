# Public release security audit

Scope: VERSION_C_CODE_RELEASE_FINAL_CANDIDATE and its eventual ZIP.

Files audited: 25
PASS_CANDIDATE: 18
AUTHOR_REVIEW_REQUIRED: 7

No patient-level rows, clinical workbook, identity map, API key, token or local absolute user path was found in the candidate scan.
CSV/GZ/derived-configuration files are flagged for human review by extension or provenance. The fixed ECM gene lists are derived configuration and their redistribution/attribution terms must be confirmed.
The Zenodo record and the GEO repositories are referenced but not rehosted. No third-party raw dataset is included.

## File-level flags

| Path | Decision | Flags | Note |
|---|---|---|---|
| `00_CONFIG/analysis_definitions.yaml` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `00_CONFIG/fixed_gene_lists/ECM33_genes.txt` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `00_CONFIG/fixed_gene_lists/ECM64_FUNCTIONAL_PANEL_MAP.csv` | AUTHOR_REVIEW_REQUIRED | NONE | Code/data-content review flag; no patient rows were found in the candidate scan. |
| `00_CONFIG/fixed_gene_lists/ECM64_genes.txt` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `01_CLINICAL_MODELS/m12_auxiliary_mar_source.py` | AUTHOR_REVIEW_REQUIRED | clinical_row_signal | Code/data-content review flag; no patient rows were found in the candidate scan. |
| `02_BULK_DISCOVERY/gse225731_pathway_discovery.R` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `03_EXTERNAL_ASSESSMENT/gse198520_fixed_components_audit.R` | AUTHOR_REVIEW_REQUIRED | clinical_row_signal | Code/data-content review flag; no patient rows were found in the candidate scan. |
| `04_XENIUM_CELL_MODULES/xenium_annotation_and_spatial_null.py` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `04_XENIUM_CELL_MODULES/xenium_multidonor_batch.py` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `06_FIGURE_REPRODUCTION/REPlot_SOURCE_RELEASE_STATUS.md` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `06_FIGURE_REPRODUCTION/figure3_source.py` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `06_FIGURE_REPRODUCTION/write_hash_manifest.R` | AUTHOR_REVIEW_REQUIRED | clinical_row_signal | Code/data-content review flag; no patient rows were found in the candidate scan. |
| `07_ENVIRONMENT/ENVIRONMENT_STATUS.md` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `07_ENVIRONMENT/R_sessionInfo_GSE225731.txt` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `07_ENVIRONMENT/R_sessionInfo_clinical_stage7r.txt` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `08_SAFE_AGGREGATE_OUTPUTS/GSE225731_PATHWAY_DISCOVERY_ALL.csv` | AUTHOR_REVIEW_REQUIRED | NONE | Code/data-content review flag; no patient rows were found in the candidate scan. |
| `08_SAFE_AGGREGATE_OUTPUTS/README.md` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `08_SAFE_AGGREGATE_OUTPUTS/Table_S1_ECM_module_results.csv` | AUTHOR_REVIEW_REQUIRED | NONE | Code/data-content review flag; no patient rows were found in the candidate scan. |
| `09_TESTS_AND_DEMO/clinical_synthetic_example.csv` | AUTHOR_REVIEW_REQUIRED | NONE | Code/data-content review flag; no patient rows were found in the candidate scan. |
| `09_TESTS_AND_DEMO/smoke_test.py` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `10_PROVENANCE/RELEASE_SCOPE.md` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `10_PROVENANCE/RELEASE_STATUS.md` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `CITATION.cff` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `LICENSE_REVIEW_REQUIRED.md` | PASS_CANDIDATE | NONE | No flagged pattern. |
| `README.md` | PASS_CANDIDATE | NONE | No flagged pattern. |
