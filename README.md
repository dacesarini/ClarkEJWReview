# Main Clark review

**Edit main.tex.** It is the main review source, with references.bib as its bibliography. Compiled output is main.pdf.

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

The manuscript's substantive claims have not been re-audited or rewritten as part of this organizational pass. Consult ../notes/review_criticisms_2026-09-12.md for the current discussion framework.
