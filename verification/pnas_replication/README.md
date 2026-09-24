# PNAS equation 3 replication audit

Status, 15 September 2026: the computations are reproducible and tested; **Clark's exact estimation procedure has not been reproduced**. Candidate weighting methods must not be described as his code. This folder does not alter source data or the review manuscript.

## Run

Requires Python, NumPy, SciPy, and `pdftotext` on PATH. Run from this directory:

```powershell
python extract_sources.py
python extract_chart_cache.py
python build_inputs.py
python sibling_diagnostics.py
python estimate.py --suite
python cross_checks.py
python test_replication.py
```

Extraction is needed only when sources change. Workbook parsing takes longer than fitting; `inputs_audit.json` caches the aggregate audit. Source SHA-256 hashes are recorded. No person-level records are copied into the results. Existing project work is preserved.

## Modify the model or estimator

```powershell
python estimate.py --cohort 1860-1919 --weight normal_log
python estimate.py --cohort 1860-1919 --weight normal_log --iterate
python estimate.py --cohort 1860-1919 --weight normal_log --model ph_fixed
python estimate.py --cohort 1780-1859 --weight N --scale level --model ph_free
python estimate.py --cohort 1780-1859 --weight equal --model genetic
python estimate.py --csv custom.csv --weight supplied --model unrestricted
```

Custom CSV required columns: `rho,N,n,d_lin`; optional `se_binary,weight`. `--omit 6` excludes zero-based row 6. No rows are omitted by default. For every fit, outputs include all coefficients, predictions, normalized weights, objective, weighted R-squared, structural bounds, and convergence. The joint sampling covariance is not estimated. These fits are descriptive, not valid pedigree-adjusted inference.

## Definitions

- Equation 3: log(rho) = a + beta*n + c*d_lin, where beta = log(b).
- `unrestricted`: three freely fitted coefficients. The reported h2=exp(a) and r=(1+m)exp(c)-1 are the theta=1 interpretation, not proof of an admissible phenotypic model.
- `ph_fixed`: theta=1, m=h2*r, with h2 and r bounded between zero and one.
- `ph_free`: theta,h2,r bounded between zero and one, m=h2*r. Here exp(a)=theta*h2, not h2. Comparisons with Clark's reported h2 use exp(a).
- `genetic`: theta=1, m and h2 estimated, c=0, r=h2*m. If theta is freed, only theta*h2 and m are identified; this engine does not pretend to separately estimate theta and h2 in that model.
- `--scale log` minimizes squared log residuals; `--scale level` minimizes squared correlation residuals. Changing scale while holding weights fixed changes the estimator.

All models retain the same n and d_lin convention: n = 1,1,2,2,3,4,5,6,7,8,9 and d_lin = 0,1,0,1,0,0,0,0,0,0,0. No average-parent/grandparent or spouse moments are included.

## Weights

Weights are normalized to mean one (which does not change estimates).

| CLI name | Unnormalized weight |
|---|---|
| equal | 1 |
| N | N |
| normal_corr | (N-3)/(1-rho^2)^2 |
| normal_log | (N-3)*rho^2/(1-rho^2)^2 |
| student_corr | (N-2)/(1-rho^2) |
| student_log | (N-2)*rho^2/(1-rho^2) |
| binary_corr | 1/se_binary^2 |
| binary_log | rho^2/se_binary^2 |
| inverse_se_corr | sqrt(normal_corr) |
| inverse_se_log | sqrt(normal_log) |
| N_rho2 | N*rho^2 |
| supplied | user-provided positive weights |

The normal approximations are reference candidates, not exact binary-outcome variances. Using N-3 here does not transform the regression's dependent variable to Fisher z. The student and inverse-SE recipes are diagnostics, not recommended inverse-variance estimators for binary correlations. Binary SEs use the empirical Pearson-correlation influence function under independent pairs; tests compare them with a multinomial delta-method calculation. They are refused for printed moments whose underlying binary margins have not been established.

`--iterate` holds weights fixed in each fitting step, then updates their rho argument using fitted correlations until the predictions change by less than 1e-11 (maximum 200 steps). This is IRLS, not direct minimization of a parameter-dependent weighted objective. Nonconverged candidates remain explicitly marked and must not be treated as solutions. Binary margins are not iteratively guessed. Goldberger's discussion concerns iterated correlation-level fitting; an iterated log fit is a separate candidate, not an attributed Clark procedure.

## Artifacts

- `sources/manifest.json`: 20 original-to-extract mappings with checksums. PDFs have page markers; DOCX equations are linearized and must be checked against originals for mathematical layout.
- `chart_cache_audit.json`: sparse indexed chart values, formulas, and external workbook links from the archived Word documents. No links are followed. The referenced author workbook `Family Tree Correlations - pnas.xlsx` is not present in the inspected project.
- `inputs_audit.json`: printed education inputs/targets, literal workbook groups, independent birth-rule groups, duplicate/overlap checks, sibling reconstruction, source checksums.
- `sibling_diagnostics.json`: alternative parent, lineage, cohort, and orientation definitions. None is silently accepted as the published sibling sample.
- `fit_results.json`, `fit_summary.csv`: specification suite (published versus literal raw samples clearly distinguished).
- `cross_checks.json`: nine-outcome checks, rounding sensitivity, leave-one-out diagnostics. Rounding ranges are not confidence intervals or rigorous worst-case bounds.
- `test_replication.py`: regression benchmarks, constraint/nesting checks, global optimizer comparison, synthetic recovery, binary SE verification, weight-scale invariance, invalid-input checks.

## Remaining replication gap

The paper and author response specify WLS but not the SE formula, exact weight vector, or software commands. Neither the published package nor the inspected referee materials supply an executable author estimation script. Archived R files are a separate correlation-level fitting exercise, not evidence of Clark's implementation. The archived `pnasrevisions.zip` has no entries.

The education workbook also has an apparent sibling/second-cousin export or labeling problem and inconsistent cohort flags. Preserve these findings; do not replace the missing sibling moments with second cousins. An exact replication needs the correct analysis sample and the original moment SEs/weights or author code. See the project-level audit note for details and the precise follow-up request.
