from pathlib import Path
import csv,json
R=Path('EJW_review');S=Path('analysis/restriction_bias_2026-09-24')
def write(n,s):(R/n).write_text(s.strip()+'\n',encoding='utf-8')
f=R/'sections/environmental_evidence.tex';s=f.read_text(encoding='utf-8').replace('FrisellEtAl2012,','');f.write_text(s,encoding='utf-8')
rows=list(csv.DictReader((S/'misspecification_population.csv').open()))
table=[]
for i,name in enumerate(dict.fromkeys(x['case'] for x in rows),1):
 rr=[x for x in rows if x['case']==name];r=rr[0]
 table.append(f"{i} & Truth & {float(r['truth_h2']):.4f} & {float(r['truth_m']):.4f} & {float(r['truth_r']):.4f} & --- \\\\")
 for r in rr:table.append(f"{i} & {'Restricted' if r['restricted']=='True' else 'Unrestricted'} & {float(r['h2']):.4f} & {float(r['m']):.4f} & {float(r['r']):.4f} & {float(r['rmse']):.4f} \\\\")
write('appendices/restriction_simulations.tex',r'''
\section{Restriction Omission and Misspecification: Controlled Experiments}
\label{app:restriction-simulations}
\subsection{Correct Specification, Identification, and Finite Samples}
If sample moments converge to $q(\psi_0)$, a minimum-distance estimator
with limiting positive-definite weight matrix $W$ has population objective
\[
 Q(\psi)=[q(\psi_0)-q(\psi)]^{\mathsf T}W[q(\psi_0)-q(\psi)].
\]
When the unrestricted moments uniquely identify $\psi_0$, this objective
has its unique minimum at the truth. Omitting an equality satisfied by
$\psi_0$ does not alter that minimum. Consistency still requires the usual
regularity conditions. Finite-sample unbiasedness does not follow.

We simulated independent Gaussian pedigrees under phenotypic and genetic
assortment, with $h^2=.60$, $m=.48$, and respectively $r=.80$ and $.288$.
The twelve fitted moments comprise eleven kinships through fourth cousins
plus spouses. Exact Wishart sample-scatter draws reproduce the sample
correlation distribution while retaining dependence among kinship estimates.
Each of $N=100,1{,}000,10{,}000$ uses 2,000 replications per model. Both
estimators minimize the same equal-weight squared correlation discrepancies;
free parameters are bounded between $10^{-6}$ and $1-10^{-6}$. Correlations
are not logged and no sample is dropped for a negative correlation.

At exact population moments, both estimators recover the true parameters
under each correctly specified model. Finite-sample biases can affect
several parameters and shrink with sample size. Under phenotypic assortment
at $N=100$, unrestricted means are $(.6037,.4670,.7977)$ and restricted
means $(.5976,.4784,.7996)$. At $N=10{,}000$ both are close to the truth.
Imposing a valid restriction does not uniformly reduce finite-sample bias:
under genetic assortment at $N=100$, the mean estimate of $r$ changes from
$.2871$ unrestricted to $.2796$ restricted, against truth $.288$.

A separate exact example isolates nonlinear recovery. With unbiased
Gaussian estimates of log-moment coefficients $(b_0,b_1,b_2)$, recovery
uses $\widehat h^2=e^{\widehat b_0}$,
$\widehat m=2e^{\widehat b_1}-1$, and
$\widehat r=2e^{\widehat b_1+\widehat b_2}-1$.
Positive normal variances make all three recovered parameters upward
biased by the exponential-moment formula, despite unbiased coefficients.
This is a stylized measurement experiment, not the exact sampling law of
logged Pearson correlations, and it does not show that a restriction
would eliminate the bias.

\subsection{When the Restriction Supplies Identification}
With only the unrestricted PH-form moments
$S=h^2(1+m)/2$ and $P=h^2(1+r)/2$, the population values
$(S,P)=(.444,.540)$ are consistent with both
$(h^2,m,r)=(.60,.48,.80)$ and $(.80,.11,.35)$, among infinitely many
triples. Imposing $m=h^2r$ gives the quadratic
$t^2-2.08t+.888=0$ for $t=h^2$; only the root $.60$ is an admissible
variance share. Omitting the restriction therefore loses identification
in this limited-moment example. An arbitrary selection along the
unrestricted solution set can be wrong about every parameter even at
infinite sample size. The example does not establish nonidentification
for Clark's full set of moments.

\subsection{Misspecified Models and Their Population Targets}
We next fit the same twelve moments under three generating processes.
The covariance matrices are built from independent Gaussian innovations,
not assigned by the fitted formulas. Population fits are checked from
three starts. Each case is also simulated with 2,000 replications at
$N=1{,}000$ and $100{,}000$, retaining all samples. Large-sample means
approach the population fitted values below, rather than the true genetic
parameters. The criterion and parameter bounds are unchanged.

Case~1 generates genetic assortment but fits phenotypic-assortment
formulas. Case~2 adds a socially transmitted component to a PH-generated
core, with variance weights $.20$ and $.80$, respectively. The cultural
component is independent of the core, has spouse correlation $.20$, and
is transmitted by parental averaging plus independent innovations of
variance $(1-.20)/2$. Spouses match independently on the core and the
cultural component, so matching is not PH on the final recorded phenotype.
The actual direct genetic variance share is $.80\times.60=.48$, genetic
spouse correlation remains $.48$, and phenotype spouse correlation is
$.80\times.80+.20\times.20=.68$. This case omits a social transmission
process and incorrectly restricts the final matching mechanism; it does
not isolate those two errors from one another.

Case~3 combines a GA-generated core of weight $.80$ with an independent
cultural component of weight $.20$, both with spouse component correlation
$.48$. The cultural component follows parental averaging plus innovation
variance $(1-.48)/2$. The fitted model assumes GA with no cultural
transmission. Here the genetic and cultural kinship patterns coincide.

\begin{table}[htbp]
\centering
\caption{Population fits under misspecification}
\label{tab:misspecification}
\begin{tabular}{llrrrr}
\hline
Case & Fit & $h^2$ & $m$ & $r$ & RMSE\\
\hline
'''+ '\n'.join(table)+r'''
\hline
\end{tabular}
\end{table}

The difference between these cases is substantive. In Case~1 all three
unrestricted fitted parameters are wrong; imposing the wrong model's
restriction worsens fit and does not restore their interpretation.
Case~2 illustrates a very close fit with substantial genetic misattribution:
the unrestricted correlation RMSE is approximately $.0021$, despite a
true genetic share of $.48$ and fitted share of $.6582$.

Case~3 gives exact observational equivalence. Its nonspousal moments are
\[
 \left(.48+.20\right)\left(\frac{1+.48}{2}\right)^n,
\]
and its spouse correlation is $.48\times.48+.20\times.48=.3264$.
A genetic-only model therefore fits every moment with $h^2=.68$, $m=.48$,
and $r=.3264$, satisfying $r=h^2m$. Enforcing the restriction cannot
separate two sources of resemblance deliberately constructed to have the
same moment pattern. This demonstrates a limit to identification, not a
claim that cultural transmission necessarily follows this particular law.

\subsection{Reproducibility and Scope}
The scripts and complete results are saved in the project's
\path{analysis/restriction_bias_2026-09-24/} directory. The files
\path{simulate.py}, \path{identification.py}, and \path{misspecification.py}
implement the three experiments. CSV files retain replication estimates,
Monte Carlo standard errors, and population targets. The report distinguishes
these correlation-level fits, including spouses, from the earlier
log-correlation experiment without spouse moments. None is a replication
of Clark's estimator or a numerical assessment of bias in his reported
estimates.
''')
# Assemble clean main manuscript. Original file is preserved in editorial snapshot.
preamble=(R/'editorial/archive_2026-09-24_before_restructure/main.tex').read_text(encoding='utf-8').split(r'\begin{document}')[0]
preamble=preamble.replace(r'\renewenvironment{itemize}{\begin{enumerate}}{\end{enumerate}}','')
a=preamble.index(r'\title{');b=preamble.index(r'\author{',a)
preamble=preamble[:a]+r'''\usepackage[hidelinks]{hyperref}
\emergencystretch=2em
\title{Review of Gregory Clark:\\ \emph{To Have and Have Not:}\\
\emph{The Determinants of Social Status in England, 1600--2026}}
'''+preamble[b:]
body='\n'+r'\begin{document}'+'\n'+r'\maketitle'+'\n'
for name in ['introduction','model_overview','measurement_identification','identification_limits','environmental_evidence','genomic_evidence','historical_application','interpretation','conclusion']:
 body+=f'\\input{{sections/{name}}}\n'
body+=r'''
\clearpage
\bibliographystyle{apalike}
\bibliography{references}
\clearpage
\appendix
\renewcommand{\thesection}{\Alph{section}}
'''
for name in ['fisher_derivations','measurement_derivations','binary_calibration','assortment_fits','restriction_simulations','quotation_dossier']:
 body+=f'\\input{{appendices/{name}}}\n'
body+=r'\end{document}'+'\n'
write('main.tex',preamble+body)
# A separate editable holding document, not part of the argumentative review.
write('editorial_holding.tex',r'''
\documentclass[12pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,natbib,enumitem}
\usepackage[T1]{fontenc}\usepackage{mathpazo}
\usepackage{xr}\externaldocument{_build/main}
\usepackage[hidelinks]{hyperref}
\emergencystretch=3em
\title{Clark Review: Editorial Holding Area}
\author{}\date{}
\begin{document}\maketitle
This companion preserves the original working list. Some items are now
integrated, some require source checks, and others remain alternatives.
The current destination and status map is in
\path{editorial/STRUCTURE_AUDIT.md}; the preserved pre-reorganization draft
is in \path{editorial/archive_2026-09-24_before_restructure/}.
These notes are not assertions that every proposed criticism has been
established. In particular, old instructions to perform work may have
been superseded by the completed analyses documented in the audit.
\input{editorial/working_notes_original}
\bibliographystyle{apalike}\bibliography{references}
\end{document}
''')
print('Manuscript and companion assembled.')
