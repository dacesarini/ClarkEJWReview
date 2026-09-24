"""Free-theta phenotypic-assortment fits to all PNAS Table 2 outcomes.

Unweighted least squares in correlation levels; 11 biological-relative
moments per outcome, no spouse observation. Source values are rounded.
"""
import json
from pathlib import Path
import numpy as np
from fit_assortment_models import DATA, RELATIVES, fit

# Columns follow the paper's Table 2; rows follow RELATIVES.
OUTCOMES = [
    ('Modern status', '1910-1997'),
    ('Log house value', '1910-1997'),
    ('Deprivation index (IMD)', '1910-1997'),
    ('Company director', '1910-1997'),
    ('Occupational status', '1780-1859'),
    ('Occupational status', '1860-1919'),
    ('Higher education', '1780-1859'),
    ('Higher education', '1860-1919'),
    ('Literacy', '1725-1869'),
]
CORRELATIONS = np.array([
    [.375,.336,.269,.168,.578,.522,.479,.326,.431],
    [.396,.360,.330,.132,.595,.512,.538,.374,.336],
    [.256,.246,.172,.057,.502,.380,.382,.228,.273],
    [.313,.278,.245,.083,.451,.362,.381,.249,.196],
    [.208,.211,.144,.061,.431,.299,.291,.173,.234],
    [.133,.147,.098,.016,.341,.274,.234,.162,.191],
    [.133,.141,.073,.069,.266,.234,.176,.186,.200],
    [.084,.084,.054,.028,.189,.203,.032,.125,.222],
    [.097,.100,.052,.053,.139,.174,.071,.085,.150],
    [.059,.069,.013,.016,.070,.111,.094,.032,.052],
    [.074,.078,.028,.021,.071,.066,.080,.014,.107],
])


def main():
    results = []
    for column, (outcome, cohort) in enumerate(OUTCOMES):
        key = outcome + ' ' + cohort
        DATA[key] = {'correlations': CORRELATIONS[:,column].tolist()}
        result = fit(key, 'phenotypic', free_theta=True)
        result.update(outcome=outcome, cohort=cohort)
        result['boundaries'] = [name for name in ('theta','h2','r_true')
                                if min(result[name], 1-result[name]) < 1e-7]
        results.append(result)
        print(key, {k:round(result[k],6) for k in
                    ('h2','m','theta','r_true','theta_h2','r_observed','raw_sse')},
              'boundary:',result['boundaries'])
    previous = json.loads(Path(__file__).with_name('assortment_fit_results.json').read_text(encoding='utf-8'))
    baseline = next(f for f in previous['fits'] if f['free_theta'])
    for key in ('h2','m','theta','r_true','raw_sse'):
        np.testing.assert_allclose(results[7][key],baseline[key],rtol=1e-6,atol=1e-8)
    report = dict(source='Clark (2023), PNAS Table 2, p. 2',
                  method='Bounded free-theta phenotypic assortment, unweighted least squares in correlation levels; multistart optimization independently checked by differential evolution.',
                  relatives=RELATIVES, fits=results,
                  caveats='No spouse moments or sampling uncertainty included. Rounded published correlations. Equal-sex and common classical-error assumptions. Boundary values are constrained optima, not evidence of certainty.')
    output = Path(__file__).with_name('all_phenotypes_free_theta_results.json')
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Saved',output)


if __name__=='__main__':
    main()
