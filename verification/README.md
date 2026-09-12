# Fisher-model section: source checks and verification

Editable source: ../sections/fisher_model.tex. The main review includes it directly.

Run python verification/verify_fisher_model.py from the repository root (or run the script by absolute path). Dependencies: Python, NumPy, SymPy. The report is written to _build/fisher_model_verification.json. These tools are not required to compile the review.

## What was checked

46 exact symbolic checks, reducing polynomials under e^2+h^2=1, include:
- Spousal covariance, preservation of phenotype variance, positive-definite spouse determinant, and Gaussian conditional cross-covariances under each PH statement.
- Marginal component covariances and phenotype variances of all eight pedigree members.
- All three spouse blocks, segregation variance, and parent/child, sibling, grandparent, avuncular, and cousin component blocks.
- All six phenotype correlations, the parent/sibling difference, and random-mating limits.
- Sensitivity of the cousin moment to additional covariance between the incoming mates' genetic deviations.

The family covariance matrix is assembled from transmission and conditional mate draws, then compared with the independently stated closed-form predictions. This checks more than substitution into the displayed formulas.

Nine simulations use h^2 in {0.1, 0.36, 0.8} and rho in {0, 0.35, 0.8}, with 250,000 independent pedigrees per setting: 2.25 million pedigrees in total. They check all 136 distinct entries of the 16-by-16 component second-moment matrix in each setting (1,224 checks), plus six phenotype moments per setting (54 checks).

All passed the prespecified six-standard-error diagnostic. Maximum absolute standardized discrepancies: 3.401649 for component moments and 2.090553 for phenotype moments. Standard errors use the exact Gaussian product-variance formula, with known zero means. These correlated diagnostics are not reported as independent hypothesis tests.

## Explicit Gaussian construction

Start one founder with X ~ N(0,I_2). Given an existing person's X, draw the mate as

X_mate = C X + R Z,

where C=rho q q', Z~N(0,I_2) is fresh, and R=I_2+(sqrt(1-rho^2)-1)q q'. Thus R R'=I_2-C C'. Each mate has marginal covariance I_2, and the spouse block is C. Their jointly Gaussian distribution has zero cross-component conditional covariance given either phenotype, establishing PH. Every person also has linear own-phenotype conditional means.

At each birth draw S~N(0,(1-rho h^2)/2) and E~N(0,1) independently of the existing family and other fresh innovations. Set A_child=(A_parent1+A_parent2)/2+S. This preserves the normal marginal component law. Draw two full siblings, independent outside mates conditional on the existing family, and a child of each pair. The eight people are two founders, their two children, the two incoming mates, and the two first cousins.

This proves existence of the specified Gaussian moment construction and provides its simulation. It is not an allele-level simulation, a proof of convergence to population-genetic equilibrium, or evidence that these restrictions describe humans.

## Source crosswalk

Sources are the checked local PDFs in ../../references/, not reconstructed quotations.

- Goldberger (1978), printed pp. 8–11: classical components, mating/variance equilibrium, Table 1, and the environmental exclusion with the spousal exception.
- Goldberger, printed p. 19: the mistake of freeing A while keeping formulas requiring A=c1*c2*m, and the need to rederive moments under alternative mating specifications.
- Goldberger, Appendix pp. A-3–A-5: additive transmission and nonzero covariance of parental environmental components with the child's additive genotype.
- Otto et al., supplied 28 July 1994 manuscript, equations (22)–(30): mate covariance block and Gaussian conditional expectations. Equation (46): parent–offspring correlation; equation (55): full siblings. Remove the cultural component and set shared-environment covariance to zero.
- Clark book draft, printed p. 75: explicitly acknowledges that r, m, h^2 are linked.
- Clark, printed pp. 296–297, Table A1 and equations A1–A3: phenotype-matching kinship formulas. The new section covers individual parent, full sibling, individual grandparent, aunt/uncle, and first cousin; it does not use the average-parent or average-grandparent rows.

Notation: the new section's rho is Clark's r and Goldberger/Otto's marital m. Its rho_A is Clark's m and Goldberger's A. A_i denotes standardized individual additive value; h A_i is the genetic contribution on the phenotype scale. q contains common parameters; X_i contains random components.

## Scope of manuscript changes

Replaced the older Fisher/Goldberger section, preserving its prior source and PDF in the dated archive. The new section corrects the old blanket claim that every non-spousal genetic–environmental covariance is zero, separates conditional linearity from PH, and separates the spouse correlation restriction from stationary segregation variance.

Dependent cross-references were updated. The later parent/sibling comparison now distinguishes an observed sign discrepancy from a formal rejection and states the common-attenuation requirement. The assertion that Clark never recognizes the constraint was removed from the replaced section. Other empirical and policy sections remain working drafts and have not received a complete source audit in this task.

## Revised criticism of the assortment discussion

The concluding discussion now distinguishes acknowledging the parameter relationship (Clark p. 75) from treating it consistently as a binding implication of the fitted model. Clark p. 16 appeals to matching on broader characteristics; pp. 74–75 appeal to measurement error in an underlying phenotype; p. 286 gives the strong-assortment genetic interpretation. The new text explains why latent measurement error preserves the restriction at the latent level, whereas a different matching rule requires a corresponding formal specification and rederived moments. It criticizes the gap between the fitted pattern and its genetic interpretation without claiming that the numerical implementation has been shown to violate the constraint. Equations and verification code are unchanged.

## Measurement-error extension

Editable source: ../sections/measurement_error.tex. It is included immediately after the Fisher-model section. The baseline section now calls its latent environmental loading e_P; the observed model uses Y=aA+eE+uU. The baseline Python verifier continues to call the latent loading e internally.

Run python verification/verify_measurement_error.py. Dependencies are the same NumPy/SymPy environment as the baseline script. The full report is written to _build/measurement_error_verification.json; measurement_error_results.json is the saved report for this revision. The baseline verifier was made importable without running simulations automatically; its command-line report was rerun and confirmed byte-for-byte identical to the previous saved result.

### Definitions and source checks

- Theta is written eta: the noise share; theta=1-eta is Clark's signal-retention/attenuation factor.
- Clark explicitly calls his theta an attenuation factor on printed p. 42. In Table A1 (printed p. 296), theta multiplies individual-relative correlations. The notation is aligned directly: theta is the attenuation factor and eta is the noise share.
- The displayed model has Var(Y)=1, a^2+e^2+u^2=1, latent h^2=a^2/theta, and rho_A=h^2*rho.
- Clark's Table A1 was visually checked against the original PDF. Its average-parent row says theta*h^2 under a table labeled correlations. A literal correlation with the arithmetic mean of the observed parents is instead theta*h^2*(1+rho)/sqrt(2*(1+theta*rho)). The text distinguishes correlations, covariances, and regression slopes and does not assume the table's average rows are interchangeable with individual-relative moments.
- Equations A1–A3 on printed p. 297 supply the logarithmic regression structure. Its intercept identifies theta*h^2. Separating noise and latent heritability additionally requires the PH restriction and enough information to identify the lineal contrast, or other suitable information. The draft derives the inversion under positive assortment. It does not claim identification from the collateral-only regression.

### Analytical and simulation checks

71 new symbolic checks passed, alongside the 46 baseline symbolic checks imported by the script. They cover the error cross-terms, observed-versus-latent scaling, spouse matrix, every table row, distant-cousin recurrence, population regression inversion, midparent correlation and slope, and nonidentification from collateral moments alone.

The Jacobian determinant of (alpha,beta,gamma) with respect to (h^2,rho,eta) is
-rho / [(1-eta)(1+rho)(1+h^2*rho)].
It is nonzero in the stated positive-assortment interior and vanishes at random mating. An exact pair of different admissible parameter triples is checked to have the same entire set of collateral moments but different parent–child moments.

Thirteen core simulations use 200,000 independently generated pedigrees each: h^2 in {0.2,0.6}, rho in {0.25,0.75}, and eta in {0,0.25,0.6}, plus (h^2,rho,eta)=(0.36,0,0.4). Each pedigree is extended through fourth cousins with fresh outside mates and births. Total: 2.6 million core pedigrees.

All 156 individual-relative moment checks, 26 midparent checks, 468 observed covariance checks, and 1,664 error-versus-true-component cross-covariance checks passed the six-standard-error diagnostic. Largest absolute standardized discrepancies were 2.515, 1.881, 2.713, and 3.545 respectively. Individual product-moment standard errors are exact Gaussian formulas; the sample midparent correlation and slope use asymptotic Gaussian standard errors. These diagnostics are not empirical specification tests.

An additional 200,000 pedigrees were used for two error counterexamples:
- Unequal reliability obeys the square-root-product attenuation rule; maximum absolute standardized discrepancy 2.061.
- Correlated errors require an additive covariance term; maximum discrepancy 2.513 for the correct formula. Applying the incorrect common-multiplier formula instead produced a discrepancy exceeding 45 standard errors.

Regression inversion was tested on exact population moments by fitting the log-linear representation and recovering the structural parameters to numerical precision. It is not a claim of finite-sample unbiasedness, and no simulated negative correlations were discarded to make a log regression appear well behaved.

### Substantive restrictions

The common factor assumes errors uncorrelated with all relevant people's true components and errors, equal signal shares in the compared samples, a common latent phenotype, and mating on that latent phenotype rather than on recording noise. The inherited restrictions excluding shared environmental transmission remain in place. Different error correlations, cohort/sex-specific reliability, or averaging relatives require rederived moments and cannot be repaired by applying a universal factor without justification.
