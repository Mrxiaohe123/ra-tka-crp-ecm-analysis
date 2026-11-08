# Reproducibility verification report

## Tests actually run

- Candidate Python syntax compilation: PASS.
- Candidate R parse check for the selected public R scripts: PASS.
- Synthetic clinical smoke test: PASS_SYNTHETIC_ONLY; this does not reproduce a clinical result.
- ZIP integrity check: PASS for the final candidate archive.

## Tests not run

The package stage did not refit the clinical LMM, CRP–VAS models, M12 MAR/MNAR
models, GSE225731 pathway associations, GSE198520 associations, Xenium donor
summaries or spatial permutations. These were intentionally not rerun under
the current scope. Public input data are not bundled, and the clinical workbook
is explicitly not shareable under the author's decision.

## Reproducibility interpretation

Syntax validity is not scientific reproduction. The public-data scripts are
conditional candidates requiring the original public inputs, source-specific
metadata and compatible package versions. The clinical scripts cannot be used
by external researchers to reproduce patient-level results because the source
workbook will not be released.

## Release gate

The candidate is not yet publishable until the author approves the code scope,
license, third-party attribution and repository location. No repository or DOI
has been created.
