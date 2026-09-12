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
