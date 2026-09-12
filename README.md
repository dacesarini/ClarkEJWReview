# Main Clark review

**Edit main.tex.** It is the main review source, with references.bib as its bibliography. Compiled output is main.pdf. The Fisher-model section is included from sections/fisher_model.tex; edit that file to revise the mathematical section.

The Git repository and its existing GitHub/Overleaf remotes remain in this directory. No remote changes were made during the folder cleanup.

## Build

Windows PowerShell: ./build.ps1

Portable sequence, from this directory:

    pdflatex main.tex
    bibtex main
    pdflatex main.tex
    pdflatex main.tex

Overleaf: main document main.tex; compiler pdfLaTeX. The bibliography uses natbib and apalike. All manuscript build dependencies are contained in this directory.

supporting/ holds the previously tracked Goldberger transcription and derivation appendix. These are standalone sources; main.tex does not automatically include them. The current, more developed Otto derivation is ../derivations/Otto_additive_unique_phenotypic_homogamy.tex.

The source-PDF catalog and full bibliography audit are in ../references/. A concise bibliography audit summary is kept here as BIBLIOGRAPHY_NOTES.md for future Git/Overleaf users.

The opening includes the eight-point criticism framework. The Fisher-model section has been replaced with an explicit, analytically and computationally checked derivation. Other sections remain working drafts; the consistency edits made alongside this replacement are not a complete audit of their substantive claims. Verification instructions and source crosswalk: verification/README.md.
