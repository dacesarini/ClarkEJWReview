# Clark review: active manuscript

The active manuscript is `main.tex`, producing [main.pdf](main.pdf). The structural audit reorganized the argument into nine main sections, with detailed derivations and evidence in six appendices. Main prose now lives in the nine files explicitly included by `main.tex`; edit those files rather than the retained legacy section sources.

- [Current section/subsection outline](editorial/OUTLINE.md)
- [Structural audit and complete to-do disposition map](editorial/STRUCTURE_AUDIT.md)
- [Editorial holding PDF](editorial_holding.pdf): all 31 original working-list items, with superseded instructions identified in the audit
- [Holding index and preserved drafts](editorial/HOLDING_INDEX.md)
- [Pre-reorganization source snapshot](editorial/archive_2026-09-24_before_restructure/main.tex)

## Build

From this folder in PowerShell:

```powershell
./build.ps1
./build.ps1 -Document editorial_holding
```

Build the main manuscript first: the holding document imports its labels. The build runs pdfLaTeX, BibTeX, and two further pdfLaTeX passes, retaining intermediate files in `_build`. `main.pdf` contains the review and scholarly appendices; the editorial to-do material is separate. For Overleaf, select `main.tex`, pdfLaTeX, natbib/apalike.

## Active files

`sections/introduction.tex`, `model_overview.tex`, `measurement_identification.tex`, `identification_limits.tex`, `environmental_evidence.tex`, `genomic_evidence.tex`, `historical_application.tex`, `interpretation.tex`, and `conclusion.tex` form the main argument.

`appendices/` contains the Fisher derivation, measurement derivation, binary calibration, alternative-mating fits, restriction/misspecification simulations, and balanced quotation dossier. The old `fisher_model.tex`, `measurement_error.tex`, `assortment_comparison.tex`, and `jencks_goldberger.tex` under `sections/` are retained legacy sources and are not included.

## Research and verification

The latest simulations and methods are in `../analysis/restriction_bias_2026-09-24/REPORT.md` and `MISSPECIFICATION.md`. Code and full replication results are beside those reports. These estimator experiments do not estimate bias in Clark's actual parameters. Earlier derivation verification remains in `verification/README.md`.

The bibliography source audit remains in `../references/metadata/audit_2026-09-23/`; unresolved source-version and missing-original issues are recorded there and in `BIBLIOGRAPHY_NOTES.md`. The structural audit did not silently fill those gaps.

The structural-edit scripts in `editorial/` record how the snapshot was transformed. Do not rerun them after editing the new draft: they regenerate text from the older snapshot.
