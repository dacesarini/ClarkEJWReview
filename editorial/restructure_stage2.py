from pathlib import Path
R=Path('EJW_review');A=R/'editorial/archive_2026-09-24_before_restructure';old=(A/'main.tex').read_text(encoding='utf-8');jg=(A/'sections/jencks_goldberger.tex').read_text(encoding='utf-8')
def cut(s,a,b):return s[s.index(a):s.index(b,s.index(a))]
def write(n,s):(R/n).write_text(s.strip()+'\n',encoding='utf-8')
conf=cut(old,r'\subsection{Confounding Is Not Limited',r'\subsection{Confounding Within')
write('sections/identification_limits.tex',r'''
\section{Good Fit, Model Restrictions, and Causal Identification}
\label{sec:identification-limits}
The preceding problems concern the connection between the model and the
recorded correlations. A further problem would remain even if that connection
were correct: matching the moments does not uniquely establish their cause.
Correct specification and identification can justify consistency under
appropriate estimation conditions; they do not guarantee finite-sample
unbiasedness. Conversely, misspecified equations can fit closely while
assigning the wrong structural meaning to their parameters.

\subsection{The Lesson of the 1970s IQ Debate}
\label{sec:1970s-precedent}
The disagreement is not new. Models fitted to overlapping IQ kinship data
in the 1970s yielded sharply different heritability estimates despite
apparently satisfactory fit. \citet{Loehlin1978} traced important differences
to assumptions about dominance, assortment, and twin environments. His
constructive recommendation was to identify the consequential assumptions
and seek external evidence about them.

\citet{Goldberger1978} distinguished problems that are easily conflated.
The Honolulu model could not uniquely identify its parameters even from
population correlations; fixing a parameter could restore an estimable
version without establishing that the chosen value was correct. Separately,
misspecified models could produce wrong structural estimates despite good
fit. His simple twin example shows how excess environmental resemblance
among MZ twins can be attributed to heritability while the fitted equations
reproduce both twin correlations exactly \citep[p.~72]{Goldberger1978}.
The issue is not merely imprecision. It is what the observations establish
about the decomposition imposed on them.

\subsection{What Omitting a Model-Implied Restriction Does}
On p.~19, Goldberger criticizes freeing spouse genetic resemblance while
retaining formulas derived under Fisher's mating scheme. Alternative
mating mechanisms may remove that restriction, but require a new derivation.
That is a warning about structural specification, not a theorem that an
unrestricted estimator must be inconsistent whenever a true equality is
left unimposed. If the unrestricted moments remain valid and uniquely
identify the true parameters, the true point still fits the population
moments exactly. Omission alone then does not change the probability limit.

Failing to impose a model-implied restriction can weaken or destroy
identification and can aggravate finite-sample bias. But persistent bias
under otherwise valid identification requires an additional argument,
such as incorrect moment equations or an invalid structural interpretation.
Goldberger's warning makes that investigation necessary; it does not
substitute for it. A statistically significant failure of a necessary
restriction challenges the joint maintained assumptions. It does not
identify which assumption failed or show that every estimate is biased.

Appendix~\ref{app:restriction-simulations} separates these cases using
population calculations and simulated pedigrees. Under the correct mating
model with identified unrestricted moments, both estimators recover the
population truth. With only parent--child and sibling moments, omitting
the phenotypic-assortment restriction instead leaves multiple triples
$(h^2,m,r)$ observationally equivalent. Under the wrong mating equations,
all three fitted parameters can converge to wrong values. Imposing the
restriction of that wrong model does not repair it.

'''+conf+r'''
\subsection{A Perfect Fit Can Misclassify Cultural Transmission}
The distinction is especially clear in a deliberately constructed
counterexample. Let a socially transmitted characteristic follow the same
parental averaging process as an additive genetic characteristic, with
independent innovations and the same degree of assortment. Both are
Gaussian and independent of one another before matching. Their phenotypic
contributions then have the same kinship-correlation pattern. A model
that includes only genetic transmission absorbs their combined variance
into its genetic parameter.

In the example in Appendix~\ref{app:restriction-simulations}, the true
direct genetic variance share is $.48$, while a cultural component
contributes $.20$. The genetic-only model reports $h^2=.68$. It reproduces
all twelve fitted correlations exactly and satisfies its own restriction
$r=h^2m$. The construction is not offered as a realistic account of every
family process. It establishes that neither perfect fit nor satisfaction
of this restriction identifies the genetic interpretation. Independent
measurements or designs on which the mechanisms differ are required.
''')
reg=cut(old,r'\subsection{Regression Towards the Mean',r'\subsection{The 1970s')
twins=cut(old,r'\subsection{Confounding Within Monozygotic Twin Pairs}',r'\section{What Modern Genomic Evidence')
a=twins.index('If\nability differences') if 'If\nability differences' in twins else -1
if a!=-1:twins=twins[:a]+r'''These findings concern the identification of a schooling effect,
not the genetic variance of schooling or earnings. Clark is right that
within-MZ comparisons need not be causal, but the limitation is well
established. Broader work on sibling comparisons likewise examines
nonshared confounding, measurement error, and the assumptions required
for causal interpretation \citep{FrisellEtAl2012,SjolanderEtAl2012,SjolanderEtAl2022}.
The appropriate response is to assess a design's identifying assumptions,
not infer that schooling has no causal effect.

'''
twins=twins.replace('% =====================================================================','')
write('sections/environmental_evidence.tex',r'''
\section{What the Environmental Comparisons Establish}
\label{sec:environmental-evidence}
Part~II is potentially more informative than another fit to the same
kinship moments: it asks how outcomes change with family circumstances.
But a finding supports genetic transmission over an environmental account
only if the competing accounts predict different findings under the
conditions studied. A null estimate for one exposure is not a test of all
social mechanisms. The intervention, population, outcome, and uncertainty
must remain explicit.

'''+reg+r'''
\subsection{Family Circumstances, Compensation, and the Relevant Contrast}
Parental loss illustrates the distinction. Clark argues that predominantly
genetic determination predicts little harm from early parental death
\citep[p.~189]{Clark2026}. Yet he also recognizes the alternative:
``Either nurture was not of great impact, or social mechanisms existed to
replace lost parents'' \citep[p.~288]{Clark2026}. If relatives or public
institutions replace an important input, a small effect of bereavement
need not imply that the input was unimportant. The comparison identifies
the consequences of loss with the responses it induces, not the effect
of removing care without replacement.

Family size and birth order require similarly specific contrasts. A
twinning-based comparison must justify why the additional birth isolates
the relevant resource change, who enters the sample, and which families'
responses it identifies. A small effect along that margin does not show
that every aspect of family resources is irrelevant. Nor does parental
symmetry by itself distinguish genetic transmission from social mechanisms
that operate through both parents. These tests should be judged against
specified alternatives rather than a single undifferentiated cultural model.

Exogenous wealth shocks offer a different source of evidence. The Georgia
land-lottery study of \citet{BleakleyFerrie2016}, for example, studies random
wealth and outcomes across generations. Such designs belong in the
assessment because they change resources without inferring that change
from kinship resemblance. Their implications remain specific to the
resource shock, institutions, and generations observed; they do not by
themselves settle the effects of education or other family inputs.

'''+twins+r'''
\subsection{Schooling Effects and Genetic Variance Are Different Estimands}
An effect of an additional year of schooling is not a heritability estimate,
and neither is the return on an additional pound of educational spending.
Evidence on compulsory-schooling reforms concerns the students induced to
remain in school by a particular reform and the outcomes measured afterward.
Its relevance to other margins of education requires an argument. Likewise,
failure to change intergenerational rank persistence does not establish a
zero effect on earnings levels, health, or well-being.

These distinctions organize the assessment of Clark's intervention-specific
evidence. The model in Part~I does not supply a missing causal estimate in
Part~II, and a particular intervention's result does not validate all of
the exclusions in Part~I. Section~\ref{sec:interpretation} returns to the
further step from causal effects to claims about optimal spending.
''')
write('sections/historical_application.tex',r'''
\section{From Marital Assortment to Historical Explanations}
\label{sec:historical-application}
Chapter~18 moves from family resemblance to a claim about how marital
assortment changes the supply of exceptional abilities. This application
requires a separate chain of evidence: the genetic nature of the trait,
the degree and history of assortment, the change in its distribution,
and the connection from that distribution to innovation. A successful
fit in Part~I would not establish all these links.

\subsection{A Classical Variance Result and a Numerical Correction}
An increase in equilibrium genetic variance under positive assortment is
a classical result of quantitative genetics, originating in Fisher's
framework and developed in subsequent work \citep{Fisher1918,Nagylaki1982}.
Clark's contribution here is the proposed historical application, whose
premises require separate evaluation. In the chapter's simplified model,
with fixed segregation variance $S$,
\[
 V_{t+1}=\frac{1+m}{2}V_t+S,
 \qquad V(m)=\frac{2S}{1-m}.
\]
Thus $V(.6)/V(0)=2.5$, not the fivefold increase stated on p.~274
\citep{Clark2026}. The fivefold result is correct for $m=.8$, the value
used in \citet[pp.~3--4]{ClarkCummins2022}. The arithmetic correction is
separate from whether fixed segregation variance is an appropriate
approximation for the historical comparison.

\subsection{Marriage Practices Do Not Directly Measure Genetic Assortment}
Observed spouse resemblance, a latent genetic correlation, and genetic
relatedness through shared ancestors are different quantities. Comparing
them across societies requires a common trait definition and an account
of measurement, partner selection, and population history. Section~\ref{sec:genomic-evidence}
explains how family-based PGIs could supply evidence about one of these
links; an inference from marriage institutions alone does not.

The cousin-marriage calculation makes the difficulty concrete. Clark
assumes random choice within the cousin pool and substitutes the cousin
correlation into a fixed-point equation for spouse genetic resemblance
\citep[pp.~278--280]{Clark2026}. Random choice with respect to the relevant
trait is an assumption, not a consequence of marrying a cousin. Moreover,
the cousin formula derived for pedigrees with independent outside mates
cannot simply be carried over to recurrent consanguineous mating: repeated
shared ancestors change the paths of relationship. This is a requirement
for a new derivation, not a claim that a corrected equilibrium has already
been established here. The calculation does not by itself demonstrate
weaker genetic assortment outside England.

\subsection{Variance, Exceptional Outcomes, and Historical Causation}
Even a justified increase in variance would not determine the number of
innovators without a model relating measured ability to exceptional
achievement. Tail probabilities depend on distributional assumptions,
and achievement also depends on opportunities to develop and apply a
trait. A comparison of equilibrium distributions is not evidence that
historical societies occupied those equilibria or that their differences
caused industrialization. These are additional empirical claims, which
should be evaluated separately from the classical variance result.
''')
# Concise interpretation with balanced source quotations; full dossier stays in appendix.
two=cut(jg,r'\subsection{Two Distinct Critiques}',r'\subsection{``Direct')
write('sections/interpretation.tex',r'''
\section{Interpreting Genetic Estimates and Policy Claims}
\label{sec:interpretation}
The preceding sections question identification and empirical support.
A separate issue remains even if a genetic decomposition is accepted:
what does it establish about social mechanisms, opportunity, and policy?
Keeping this conditional argument separate prevents an objection to the
estimates from substituting for an objection to their interpretation.

'''+two+r'''
\subsection{Direct Genetic Effects Can Operate through Social Responses}
Clark describes the direct path as one in which ``the genome dictates
child social aptitudes, independent of any interaction with the parents''
\citep[p.~286]{Clark2026}. This excludes pathways that the modern direct-effect
estimand includes, as Section~\ref{sec:model-components} explains. An effect
of the child's own genotype can operate through parental attention,
teachers, or peers. Evidence against parental genetic nurture would not
establish that these responses are absent. Confounding asks whether an
effect has been identified; mediation asks how that effect operates.
Treating ``direct'' as narrowly physiological conflates the two.

\subsection{Clark's Qualifications and His Stronger Conclusions}
The book also contains explicit acknowledgments of these distinctions.
Clark writes: ``This does not in itself imply that social interventions
cannot change social outcomes'' \citep[p.~80]{Clark2026}, and recognizes
that intervention remains an empirical question even without convincing
evidence of genetic nurture \citep[p.~86]{Clark2026}. He notes that
heritability depends on the range of family environments, and allows
that replacement care could explain small effects of parental loss.
The criticism cannot fairly be that he never recognizes such possibilities.

The problem is their inconsistent application. The claim that direct
genetic transmission implies excessive education spending (p.~9), or
that social outcomes can change only when genetic endowments change
(p.~10), does not follow from a variance decomposition. Nor does a
genetic explanation establish meritocracy (p.~7): discriminatory social
responses to inherited characteristics can themselves generate genetic
associations or effects. Appendix~\ref{app:quotations} preserves ten
critical passages alongside ten qualifications, with their contexts.
These are tensions in the interpretation, not proof that every empirical
claim in the book is false.

\subsection{Education Spending: Effects, Costs, and the Policy Objective}
Clark calls excessive education expenditure an implication of direct
genetic transmission and points to persistent mobility despite narrowing
gaps in schooling \citep[p.~9]{Clark2026}. But whether expenditure is
excessive requires marginal benefits and costs. The relevant benefits can
include absolute earnings, health, and consumption as well as changes in
inequality or relative mobility. A program that raises everyone's reading
ability can be valuable even if the same families retain their relative
positions. Unchanged ranks or intergenerational correlations therefore do
not establish that an intervention was ineffective or not worth its cost.

The intervention-specific evidence discussed in
Section~\ref{sec:environmental-evidence} is relevant to this question, but
its estimands must be matched to the policy claim. A local effect of a
school-leaving reform, an average return to schooling, and the return to
an expansion of public expenditure need not coincide. The genetic model
cannot supply the missing welfare comparison.

\subsection{Historical and Population Extrapolations Need Additional Evidence}
The same logic applies to claims about future outcomes and population
composition. Genetic contributions estimated under current institutions
do not establish invariant effects under different schools, labor markets,
or tax systems. Clark's acknowledgment that social and legal restrictions
limited the prospects of Jewish migrants in their countries of origin
\citep[p.~291, n.~105]{Clark2026} illustrates the relevance of those
conditions. Conclusions about future fiscal costs or historical innovation
require their own evidence; they cannot be read directly from a fitted
kinship model. Neither the identification critique nor the interpretive
critique implies that any particular policy must succeed. They establish
why success or failure must be assessed on evidence appropriate to it.
''')
write('sections/conclusion.tex',r'''
\section{Conclusion}
The book supplies valuable genealogical data, documents persistent family
resemblance, and raises consequential questions about its causes. Its
attention to measurement is also well placed. These contributions do not
depend on accepting a single explanation of the correlations.

The genetic interpretation requires more than close fit. The mating
mechanism must generate the equations being estimated, the measurement
model must connect those equations to the recorded outcomes, and the
identifying exclusions must be confronted with evidence. Environmental
comparisons can help, but only against alternatives that predict different
results. Modern genomic designs are especially valuable because random
transmission supplies information unavailable in average kinship
correlations. They can support the additive approximation while challenging
other exclusions.

Even valid genetic estimates would not divide social and genetic causation
into mutually exclusive mechanisms or establish the effects and value of
interventions. Clark recognizes parts of this argument, but his stronger
historical and policy conclusions do not consistently respect it. The
review's disagreement is therefore with the steps connecting evidence to
causal explanation and policy, not with the possibility that genes matter.
''')
print('Argument sections drafted.')
