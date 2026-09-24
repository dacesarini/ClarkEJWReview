"""Attempt the described log regression, and audit supplied education pairs.

Not an exact replication: Clark's per-moment standard errors and executable
analysis code have not been located. Weight alternatives are explicitly named.
Raw groups are interpreted literally using their labels and cohort flags.
"""
import json
import math
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np
from fit_assortment_models import DATA, RELATIVES

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
LABELS = ['sons', 'father-son', 'siblings-rem', 'grandson', 'cousin',
          'cousin-rem', 'cousin2', 'cousin2-rem', 'cousin3', 'cousin3-rem', 'cousin4']
DISTANCE = np.array([1,1,2,2,3,4,5,6,7,8,9])
LINEAL = np.array([0,1,0,1,0,0,0,0,0,0,0])
DESIGN = np.column_stack([np.ones(11), DISTANCE, LINEAL])


def regression(y, weights):
    y, weights = np.array(y), np.array(weights, dtype=float)
    root = np.sqrt(weights / weights.mean())
    a, logb, c = np.linalg.lstsq(DESIGN*root[:, None], np.log(y)*root, rcond=None)[0]
    product = np.exp(a)
    m = 2*np.exp(logb)-1
    r = (1+m)*np.exp(c)-1
    h2 = m/r
    theta = product/h2
    return dict(a=float(a),log_b=float(logb),c=float(c),exp_a=float(product),
                m=float(m),r=float(r),recovered_h2=float(h2),recovered_theta=float(theta),
                admissible_free_theta_PH=bool(0<h2<=1 and 0<theta<=1 and 0<=r<=1 and 0<=m<=h2))


def read_cells(path):
    cells = defaultdict(lambda: [0,0,0,0])
    with zipfile.ZipFile(path) as z:
        strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            strings = [''.join(t.text or '' for t in row.findall('.//m:t',NS))
                       for row in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        # Resolve the sheet by name, not by relying on a fixed sheet number.
        workbook = ET.fromstring(z.read('xl/workbook.xml'))
        links = {r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheet = next(s for s in workbook.find('m:sheets',NS) if s.attrib['name']=='Table 2 Ded 1780-1919')
        target = links[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
        target = target.lstrip('/') if target.startswith('/') else 'xl/'+target
        with z.open(target) as stream:
            for _, row in ET.iterparse(stream, events=('end',)):
                if not row.tag.endswith('}row'):
                    continue
                values = {}
                for cell in row:
                    value = cell.find('m:v',NS)
                    if value is None:
                        continue
                    text = value.text
                    if cell.attrib.get('t')=='s':
                        text = strings[int(text)]
                    values[''.join(c for c in cell.attrib['r'] if c.isalpha())] = text
                if row.attrib.get('r')=='1':
                    if [values.get(c) for c in ('A','E','F','I','M')] != ['relationship','per1780','per1860','ded0','ded1']:
                        raise ValueError('Unexpected source columns')
                elif 'I' in values and 'M' in values:
                    x,y = float(values['I']),float(values['M'])
                    if x not in (0,1) or y not in (0,1):
                        raise ValueError('Nonbinary education outcome')
                    for cohort, column in [('1780-1859','E'),('1860-1919','F')]:
                        if values.get(column)=='1':
                            cells[(cohort,values.get('A'))][int(x)*2+int(y)] += 1
                row.clear()
    return cells


def binary_stats(cells):
    count = sum(cells)
    px,py = (cells[2]+cells[3])/count,(cells[1]+cells[3])/count
    sx,sy = math.sqrt(px*(1-px)),math.sqrt(py*(1-py))
    rho = (cells[3]/count-px*py)/(sx*sy)
    influence_variance = 0
    for index,frequency in enumerate(cells):
        x,y = index//2,index%2
        zx,zy = (x-px)/sx,(y-py)/sy
        influence = zx*zy-rho*(zx*zx+zy*zy)/2
        influence_variance += frequency/count*influence**2
    return dict(n=count,correlation=rho,binary_iid_se=math.sqrt(influence_variance/count),cells=cells)


def main():
    project = Path(__file__).resolve().parents[2]
    raw = read_cells(project/'Data/pnas.2300926120.sd03.xlsx')
    report = dict(note='Literal raw-label audit, not a repaired dataset or exact replication. Binary SEs assume independent pairs; normal approximations are alternative assumptions, not documented Clark weights.',cohorts={})
    for cohort in DATA:
        y = np.array(DATA[cohort]['correlations'])
        n = np.array(DATA[cohort]['pairs'])
        variants = {'OLS':np.ones(11),'pair_counts':n,
                    'normal_approx_inverse_var_correlation':(n-1)/(1-y*y)**2,
                    'normal_approx_inverse_var_log_correlation':(n-1)*y*y/(1-y*y)**2}
        rows = [dict(relationship=name,raw_label=label,**binary_stats(raw[(cohort,label)]),
                     published_correlation=float(obs),published_n=int(total))
                for name,label,obs,total in zip(RELATIVES,LABELS,y,n)]
        raw_y = np.array([r['correlation'] for r in rows])
        raw_se = np.array([r['binary_iid_se'] for r in rows])
        report['cohorts'][cohort] = dict(
            published_rounded_regressions={name:regression(y,w) for name,w in variants.items()},
            raw_label_audit=rows,
            raw_regressions={'inverse_binary_iid_var_correlation':regression(raw_y,1/raw_se**2),
                             'inverse_binary_iid_var_log_correlation':regression(raw_y,raw_y**2/raw_se**2)})
        print(cohort,flush=True)
        for row in rows:
            if row['n']!=row['published_n'] or round(row['correlation'],3)!=row['published_correlation']:
                print('Mismatch:',row,flush=True)
        print('Raw regression attempts:',report['cohorts'][cohort]['raw_regressions'],flush=True)
    report['book_rounded_OLS_only'] = {
        '1780-1859':regression([.558,.500,.398,.426,.366,.235,.103,.098,.186,.118,.020],np.ones(11)),
        '1860-1919':regression([.359,.305,.213,.273,.146,.181,.045,.076,.099,.079,.040],np.ones(11))}
    path = Path(__file__).with_name('appendix_regression_audit.json')
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Saved',path,flush=True)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
