param([ValidateSet('main','editorial_holding')][string]$Document = 'main')
$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    New-Item -ItemType Directory -Path '_build' -Force | Out-Null
    if ($Document -eq 'editorial_holding') {
        if (-not (Test-Path -LiteralPath '_build/main.aux')) {
            throw 'Build main.tex before the editorial holding document.'
        }
        Get-Content -LiteralPath '_build/main.aux' | Where-Object { $_.StartsWith('\newlabel') } | Set-Content -LiteralPath '_build/main-labels.aux' -Encoding UTF8
    }
    function Invoke-Checked([string]$Program, [string[]]$Arguments) {
        & $Program @Arguments
        if ($LASTEXITCODE -ne 0) { throw "$Program failed with exit code $LASTEXITCODE" }
    }
    Invoke-Checked 'pdflatex' @('-interaction=nonstopmode','-halt-on-error','-output-directory=_build',"$Document.tex")
    Invoke-Checked 'bibtex' @("_build/$Document")
    Invoke-Checked 'pdflatex' @('-interaction=nonstopmode','-halt-on-error','-output-directory=_build',"$Document.tex")
    Invoke-Checked 'pdflatex' @('-interaction=nonstopmode','-halt-on-error','-output-directory=_build',"$Document.tex")
    Copy-Item -LiteralPath "_build/$Document.pdf" -Destination "$Document.pdf" -Force
} finally { Pop-Location }
