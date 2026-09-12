$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    New-Item -ItemType Directory -Path '_build' -Force | Out-Null
    function Invoke-Checked([string]$Program, [string[]]$Arguments) {
        & $Program @Arguments
        if ($LASTEXITCODE -ne 0) { throw "$Program failed with exit code $LASTEXITCODE" }
    }
    Invoke-Checked 'pdflatex' @('-interaction=nonstopmode','-halt-on-error','-output-directory=_build','main.tex')
    Invoke-Checked 'bibtex' @('_build/main')
    Invoke-Checked 'pdflatex' @('-interaction=nonstopmode','-halt-on-error','-output-directory=_build','main.tex')
    Invoke-Checked 'pdflatex' @('-interaction=nonstopmode','-halt-on-error','-output-directory=_build','main.tex')
    Copy-Item -LiteralPath '_build/main.pdf' -Destination 'main.pdf' -Force
} finally { Pop-Location }
