from pathlib import Path
R=Path('EJW_review');A=R/'editorial/archive_2026-09-24_before_restructure'
old=(A/'main.tex').read_text(encoding='utf-8');f=(A/'sections/fisher_model.tex').read_text(encoding='utf-8');me=(A/'sections/measurement_error.tex').read_text(encoding='utf-8');jg=(A/'sections/jencks_goldberger.tex').read_text(encoding='utf-8');assort=(A/'sections/assortment_comparison.tex').read_text(encoding='utf-8')
def cut(s,a,b):return s[s.index(a):s.index(b,s.index(a))]
def write(n,s):(R/n).write_text(s.strip()+'\n',encoding='utf-8')
write('sections/introduction.tex',r'''
\section{Introduction}
Gregory Clark's \emph{To Have and Have Not} assembles an unusually rich
record of family resemblance and social mobility in England. Its Families
of England and Marriages of England databases connect outcomes across
many generations and types of relatives. These data, and the book's
attention to the measurement of social status, are substantial contributions
independent of the explanation Clark favors \citep{Clark2026}.

The book's central claim is more demanding than persistence itself.
Clark interprets the pattern of kinship correlations as evidence that
additive genetic transmission, combined with strong assortative mating,
accounts for social resemblance across generations. Part~II examines
family and environmental comparisons that he takes to distinguish this
account from social transmission. Part~III extends the argument to
education, historical development, and population composition.

This review asks whether those inferences follow. Three questions must
be separated. First, do the fitted equations consistently represent the
model Clark invokes, including its mating and measurement assumptions?
Second, do the data distinguish that model from alternative processes
that generate family resemblance? Third, even if the genetic estimates
were valid, what would they establish about social mechanisms and the
benefits of intervention? A defect at one stage does not establish every
criticism at the others, and accepting genetic influences does not settle
any of these questions.

Measurement deserves particular care. Noisy occupational and educational
proxies can distort comparisons of mobility, but labeling the difference
between a recorded outcome and a latent ability ``measurement error''
does not establish the classical error assumptions required by a correction.
Nor can the consequences of recording a binary outcome generally be
represented by one common attenuation factor. These issues affect what
the reported parameters measure before their causal interpretation arises.

The argument proceeds from the model's assumptions to measurement and
estimation, then to identification and independent evidence. The discussion
of modern genomics explains how random genetic transmission can test
assumptions previously imposed on kinship correlations. Separate sections
assess the historical application of assortment and the interpretation of
genetic estimates for policy. Detailed derivations, simulations, and a
balanced dossier of quotations are retained in the appendices. The central
criticism is not that the observed persistence is uninteresting or that
genetic explanations are inadmissible. It is that the evidence does not
establish the particular genetic explanation, or its wider implications,
as securely as the book claims.
''')
env=cut(f,r'\subsection{The Strength of the Environmental Covariance Exclusions}',r'\subsection{What Phenotypic Homogamy Assumes}')
env=env.replace('In the birth\nprocess below','In the birth\nprocess of Appendix~\\ref{app:fisher}').replace('the derivation below gives','the derivation in Appendix~\\ref{app:fisher} gives').replace('with siblings\' genes','with siblings\' genes')
env=env.replace(r'$\operatorname{Cov}(E_M,A_o)=\rho h e_P/2$',r'$\operatorname{Cov}(E_M,A_o)=r h e_P/2$')
write('sections/model_overview.tex',r'''
\section{What the Genetic Model Assumes and Predicts}
\label{sec:fisher-model}
\subsection{The Outcome and the Genetic Estimand}
\label{sec:model-components}
The benchmark is an additive model for a specified phenotype, not a general
theory of every outcome called social status. Write the standardized phenotype
as
\[
 P_i=hA_i+e_P E_i,\qquad e_P=\sqrt{1-h^2},
\]
where $A_i$ and $E_i$ have zero means, unit variances, and zero within-person
covariance. The additive genetic contribution is $hA_i$, with variance
$h^2$. Define $m=\operatorname{Corr}(A_M,A_F)$ and
$r=\operatorname{Corr}(P_M,P_F)$ for spouses. The sexes share component
variances, and those variances are stationary across generations.
Appendix~\ref{app:fisher} gives the primitive assumptions and derivations;
there $\rho=r$ and $\rho_A=m$.

The intended causal genetic contribution must be distinguished from a
population predictor. In modern genomics, a direct effect originates in
the individual's own genotype; an indirect effect originates in another
person's genotype, apart from transmission of the individual's alleles.
The direct effect can include parental, teacher, or peer responses to the
individual's genetically influenced characteristics. It is not an effect
that holds all social responses fixed
\citep[Sections~3.1 and~3.4]{BenjaminEtAl2026}. This definition matters both
for the covariance exclusions and for Clark's later interpretation of
``social'' effects.

'''+env+r'''
\subsection{The Mating Mechanism Restricts the Parameters}
Under primary phenotypic assortment, spouses match on $P$ rather than
separately on its genetic and environmental components. With the conditional
linearity and independence assumptions stated in Appendix~\ref{app:fisher},
\begin{equation}
 m=h^2r.\label{eq:overview-restriction}
\end{equation}
The relation is derived from the matching mechanism; it is not an optional
additional assumption or a consequence that requires an intergenerational
equilibrium argument. A high phenotype predicts a genetic value of $hP$;
applying that prediction on both sides of a marriage gives the factor $h^2$.

Assortment on genetic value instead gives $r=h^2m$ under independent
residual environments. Assortment on several traits can generate still
other restrictions. Clark's appeal to matching on ``something deeper''
therefore requires a specification: either it is the same phenotypic
model with a measurement equation, or it changes the mating process and
requires new kinship formulas. These possibilities cannot be exchanged
while retaining the original parameter interpretation.

\subsection{Inheritance and the Predictions for Relatives}
The child receives half each parent's additive genetic value plus a
mean-zero segregation deviation:
\[
 A_o=\tfrac12(A_M+A_F)+S_o.
\]
Independent segregation across births and stationary unit genetic variance
require $\operatorname{Var}(S_o)=(1-m)/2$. The child's environmental draw
is fresh in the sense specified above. Under phenotypic assortment the
resulting correlations are
\[
 \operatorname{Corr}(P_p,P_o)=\frac{h^2(1+r)}2,
 \qquad
 \operatorname{Corr}(P_1,P_2)=\frac{h^2(1+m)}2
\]
for a parent and child and for full siblings, respectively. Thus
\begin{equation}
 \operatorname{Corr}(P_p,P_o)-\operatorname{Corr}(P_1,P_2)
 =\frac{h^2r(1-h^2)}2>0
 \label{eq:overview-gap}
\end{equation}
when assortment is positive and $0<h^2<1$. More distant moments require
additional assumptions about incoming spouses and shared ancestors,
not merely Mendelian transmission. The complete derivation exposes where
each exclusion enters. The next question is whether these latent
predictions survive the way Clark's outcomes are recorded and fitted.
''')
# Technical derivation: retain proofs, remove duplicated main-text critique and environmental essay.
fa=f.replace(r'\section{The Fisher Model with Phenotypic Assortment}',r'\section{Derivation of the Phenotypic-Assortment Model}').replace(r'\label{sec:fisher-model}',r'\label{app:fisher}').replace(r'\label{sec:model-components}',r'\label{app:model-components}')
a=fa.index(r'\subsection{The Strength of');b=fa.index(r'\subsection{What Phenotypic',a)
fa=fa[:a]+r'''The substantive interpretation of the environmental exclusions is in
Section~\ref{sec:environmental-exclusions}. We retain those exclusions
throughout this derivation.

'''+fa[b:]
a=fa.index(r'\citet[p.~19]{Goldberger1978}');b=fa.index(r'\paragraph{Verification.}',a)
fa=fa[:a]+r'''The implications for identification, specification testing, and bias
are discussed in Section~\ref{sec:identification-limits}. An unrestricted
moment fit can be used to test restrictions without constituting a fit
of the restricted structural model.

'''+fa[b:]
write('appendices/fisher_derivations.tex',fa)
# Separate binary calibration from measurement derivation.
a=me.index(r'\paragraph{A binary-outcome calibration.}');b=me.index(r'\subsection{Averaging Parents',a)
binary_details=me[a:b].replace(r'\paragraph{A binary-outcome calibration.}\mbox{}\par',r'\section{Binary-Outcome Calibration: Design and Results}\label{app:binary-calibration}',1)
write('appendices/binary_calibration.tex',binary_details)
ma=me[:a]+me[b:]
ma=ma.replace(r'\section{Adding Classical Measurement Error}',r'\section{Derivation of the Measurement-Error Correction}').replace(r'\label{sec:measurement-error}',r'\label{app:measurement-error}')
write('appendices/measurement_derivations.tex',ma)
# Main measurement/identification section, with one binary discussion.
binmain=cut(old,r'\subsection{Binary Outcomes Are Not Continuous Phenotypes}',r'\subsection{Confounding Is Not Limited')
a=binmain.index('The issue is quantitative');b=binmain.index('The continuous version',a)
binmain=binmain[:a]+r'''A calibrated Gaussian-threshold example shows the practical importance.
The latent generating parameters are $h^2=.60$, $r=.80$, $m=.48$, and
continuous reliability $\theta=.80$. We threshold the recorded continuous
outcome at the approximately 2.5\% educational prevalence in Clark's
historical cohorts. Appendix~\ref{app:binary-calibration} reports the design,
prevalence construction, and complete calibration table.

'''+binmain[b:]
a=binmain.index('Across the grid');b=binmain.index('Thresholding changes',a)
binmain=binmain[:a]+binmain[b:]
a=binmain.index('The same pattern is visible');b=binmain.index('These simulations are not',a)
binmain=binmain[:a]+r'''This problem concerns the observation equation even if the latent
inheritance model is correct. A threshold model can relate observed
prevalences and joint probabilities to latent correlations. One freely
estimated common attenuation factor cannot generally substitute for that
mapping; a separate free factor for every kinship would instead surrender
the cross-kinship restrictions being tested.

'''+binmain[b:]
empirical=cut(old,r'\subsection{The Same Pattern, on Clark',r'\subsection{Binary Outcomes')
empirical=empirical.replace(r'\subsection{The Same Pattern, on Clark\'s Own Data}',r'\subsection{A Diagnostic in Clark\'s Reported Correlations}') if False else empirical.replace('The Same Pattern, on Clark\'s Own Data','A Diagnostic in Clark\'s Reported Correlations')
empirical=empirical.replace(r'Section~\ref{sec:model-tests}',r'Section~\ref{sec:fisher-model}').replace(r'\eqref{eq:gap}',r'\eqref{eq:overview-gap}')
empirical=empirical.replace('The higher-education point estimates have the\nopposite sign in both periods.', 'The higher-education point estimates have the\nopposite sign in both periods. Their binary scale and potentially different\nkinship-sample prevalences require the additional observation model just\ndiscussed; this table alone is not a test of a latent-liability restriction.')
# Preserve the illustrative assortment fit in an appendix, summarize result in the main text.
assort=assort.replace(r'\section{Assortment, Measurement Error, and What the Fitted Parameters Identify}',r'\section{Illustrative Fits under Alternative Mating Models}').replace(r'\label{sec:assortment-comparison}',r'\label{app:assortment-comparison}')
write('appendices/assortment_fits.tex',assort)
write('sections/measurement_identification.tex',r'''
\section{From Observed Correlations to Genetic Parameters}
\label{sec:measurement-error}\label{sec:assortment-comparison}
\subsection{What Classical Measurement Error Assumes}
Let the standardized recorded outcome be
\[
 Y_i=\sqrt{\theta}P_i+\sqrt{1-\theta}U_i,
 \qquad 0<\theta\leq1.
\]
Here $\theta$ is reliability, not the noise share. If recording errors are
uncorrelated with all relevant relatives' true components and with one
another, and reliability is common, then
\[
 \operatorname{Corr}(Y_i,Y_j)=\theta\operatorname{Corr}(P_i,P_j)
 \quad(i\ne j).
\]
Attenuation occurs at the two observed endpoints, not at each intervening
generation. The genetic variance share of $Y$ is $\theta h^2$, whereas
$h^2$ refers to $P$. Appendix~\ref{app:measurement-error} derives these
relations and the different normalization required for parental averages.

These assumptions go beyond imperfect measurement. Shared classification
errors add covariance; unequal reliability produces the pair-specific
factor $\sqrt{\theta_i\theta_j}$. More fundamentally, the difference between
education and an asserted latent social ability can contain family resources,
institutional access, and preferences, rather than independent recording
noise. A plausible value of the fitted attenuation factor does not validate
these exclusions.

\subsection{Identification Depends on the Mating Model}
For the nonspousal kinships in Clark's regression, the phenotypic-assortment
form is
\[
 q_{n,D}=\theta h^2\left(\frac{1+m}{2}\right)^n
 \left(\frac{1+r}{1+m}\right)^D,
\]
where $n$ indexes genealogical distance and $D=1$ for the lineal rows.
The three coefficients of its logarithm identify $\theta h^2$, $m$, and
$r$, when the data supply the necessary distance and lineal variation.
Under phenotypic assortment, the restriction $m=h^2r$ then separates
$h^2=m/r$ from $\theta$. With known attenuation, that same restriction
instead reduces the number of freely fitted parameters. Arbitrary
regression coefficients need not yield admissible structural values.

Under genetic assortment, the lineal adjustment disappears. These
kinship correlations identify $m$ and the product $\theta h^2$, but not
reliability and heritability separately. Even the observed spouse
correlation is then $\theta h^2m$. The intuition is that untransmitted
true environmental variation and recording noise both dilute the genetic
signal without changing mate selection on genetic value. They cannot be
separated by naming one of them measurement error.

The book's equation~A3 therefore requires a precise interpretation. Three
regression coefficients can be legitimate under phenotypic assortment
with unknown attenuation, provided recovery invokes the restriction and
checks admissibility. Their number alone is not an error. But the
intercept identifies $\theta h^2$, not latent heritability, and changing
the mating mechanism changes the permissible recovery. Illustrative fits
in Appendix~\ref{app:assortment-comparison} show how similar fit can accompany
different identification claims. They are descriptive comparisons, not
formal tests using a joint sampling covariance matrix.

'''+binmain+empirical)
write('sections/genomic_evidence.tex',cut(old,r'\section{What Modern Genomic Evidence Can Identify}',r'\input{sections/jencks_goldberger}').replace('as discussed below','as defined in Section~\\ref{sec:model-components}'))
# Keep the exact full quotation dossier in one appendix, without duplicating it in the main discussion.
quotes=cut(jg,r'\subsection{Passages Vulnerable',r'\subsection{The Inconsistency')
write('appendices/quotation_dossier.tex',r'\section{Quotations and Qualifications in Clark}\label{app:quotations}'+'\n'+quotes)
# Preserve complete original working notes as a separate editorial document body.
working=cut(old,r'\section*{Working Notes: Points to Revisit}',r'\end{document}')
write('editorial/working_notes_original.tex',working)
# Preserve original sections no longer used, as documented source material in the snapshot.
print('Core sections and appendices written; original sources retained.')
