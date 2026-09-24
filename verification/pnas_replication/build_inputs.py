"""Audit education inputs. Literal data, rule-based alternatives, no silent repairs."""
from collections import defaultdict, Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from workbooks import rows, records

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
COHORTS = {'1780-1859':(1780,1860), '1860-1919':(1860,1920)}
LABELS = ['sons','father-son','siblings-rem','grandson','cousin','cousin-rem',
          'cousin2','cousin2-rem','cousin3','cousin3-rem','cousin4']
NAMES = ['Full sibling','Child','Sibling-rem','Grandchild','Cousin','Cousin-rem',
         'Cousin2','Cousin2-rem','Cousin3','Cousin3-rem','Cousin4']
DISTANCE = [1,1,2,2,3,4,5,6,7,8,9]
LINEAL = [0,1,0,1,0,0,0,0,0,0,0]
PUBLISHED = {
    '1780-1859':dict(rho=[.479,.538,.382,.381,.291,.234,.176,.032,.071,.094,.080],
                     N=[10772,8337,12030,4619,14635,12319,15935,12212,15349,12503,15953],
                     clark=dict(b=.774,m=.548,h2=.627,R2=.983)),
    '1860-1919':dict(rho=[.326,.374,.228,.249,.173,.162,.186,.125,.085,.032,.014],
                     N=[15806,15394,34874,15661,26829,51001,33813,57201,33900,53447,30104],
                     clark=dict(b=.766,m=.531,h2=.420,R2=.886))}

def binary_stats(cells):
    cells=np.asarray(cells,dtype=float);N=cells.sum()
    if N<=0: return dict(N=0,cells=cells.tolist())
    p=cells/N;x=np.array([0,0,1,1]);y=np.array([0,1,0,1])
    px=p@x;py=p@y;sx=math.sqrt(px*(1-px));sy=math.sqrt(py*(1-py))
    if sx*sy==0: return dict(N=int(N),cells=cells.astype(int).tolist(),rho=None,se_binary=None)
    zx=(x-px)/sx;zy=(y-py)/sy;rho=float(p@(zx*zy))
    influence=zx*zy-rho*(zx*zx+zy*zy)/2
    var=float(p@influence**2/N)
    return dict(N=int(N),cells=cells.astype(int).tolist(),rho=rho,px=float(px),py=float(py),se_binary=math.sqrt(var))

def pair_audit():
    path=ROOT/'Data/pnas.2300926120.sd03.xlsx'
    counts=defaultdict(lambda:np.zeros(4,dtype=int));flags=Counter();totals=Counter()
    pair_ids=defaultdict(set);duplicate=Counter();birth_missing=Counter()
    for number,row in rows(path,'Table 2 Ded 1780-1919'):
        if number==1:
            assert [row.get(c) for c in ('A','E','F','I','M')]==['relationship','per1780','per1860','ded0','ded1']
            continue
        if not ('I' in row and 'M' in row):continue
        label=row['A'];x=float(row['I']);y=float(row['M']);assert x in (0,1) and y in (0,1)
        cell=int(x)*2+int(y);totals[label]+=1
        pid=(row.get('G'),row.get('K'))
        if pid in pair_ids[label]:duplicate[label]+=1
        pair_ids[label].add(pid)
        for cohort,(lo,hi) in COHORTS.items():
            flag=row.get('E' if lo==1780 else 'F')=='1'
            if flag:counts[('flags',cohort,label)][cell]+=1
            if 'H' not in row or 'L' not in row:
                birth_missing[label]+=1;continue
            birth=(lo<=float(row['H'])<hi or lo<=float(row['L'])<hi)
            if birth:counts[('birth_rule',cohort,label)][cell]+=1
            if birth!=flag:flags[(cohort,label,'birth_only' if birth else 'flag_only')]+=1
    overlaps=[]
    for a,b in itertools.combinations(sorted(pair_ids),2):
        overlap=len(pair_ids[a]&pair_ids[b])
        if overlap:overlaps.append(dict(label_a=a,label_b=b,identical_ordered_pair_ids=overlap))
    return dict(group_stats=[dict(rule=rule,cohort=c,raw_label=l,**binary_stats(v)) for (rule,c,l),v in sorted(counts.items())],
                label_totals=dict(totals),duplicate_ordered_pairs=dict(duplicate),cross_label_overlap=overlaps,
                missing_births=dict(birth_missing),flag_disagreement=[dict(cohort=c,label=l,direction=d,N=n) for (c,l,d),n in sorted(flags.items())])

def rebuild_siblings():
    """Alternative definitions reported separately; shared father != proven full sibling."""
    groups=defaultdict(list);ids=set();duplicates=0
    for _,r in records(ROOT/'Data/pnas.2300926120.sd01.xlsx','Data'):
        if r.get('dfem')!='0' or r.get('ded') not in ('0','1') or 'byr' not in r:continue
        pid=r['pid']
        if pid in ids:duplicates+=1
        ids.add(pid)
        father=r.get('pid_fath')
        if not father or father=='0':continue
        if not 1600<=float(r['byr'])<1920:continue
        groups[father].append(r)
    counts=defaultdict(lambda:np.zeros(4,dtype=int))
    for group in groups.values():
        for a,b in itertools.combinations(sorted(group,key=lambda x:float(x['pid'])),2):
            cell=int(a['ded'])*2+int(b['ded']);same_mother=a.get('pid_moth') and a.get('pid_moth')!='0' and a.get('pid_moth')==b.get('pid_moth')
            for cohort,(lo,hi) in COHORTS.items():
                if not (lo<=float(a['byr'])<hi or lo<=float(b['byr'])<hi):continue
                counts[('shared_father',cohort)][cell]+=1
                if same_mother:counts[('same_known_parents',cohort)][cell]+=1
                if int(float(a['nid']))<1000 and int(float(b['nid']))<1000:
                    counts[('shared_father_elite',cohort)][cell]+=1
                    if same_mother:counts[('same_known_parents_elite',cohort)][cell]+=1
    return dict(duplicate_eligible_person_ids=duplicates,orientation='ascending numeric person ID, each unordered pair once',
                estimates=[dict(definition=d,cohort=c,**binary_stats(v)) for (d,c),v in sorted(counts.items())])

def main():
    print('Auditing pair workbook...',flush=True);audit=pair_audit()
    print('Rebuilding sibling candidates from person workbook...',flush=True);siblings=rebuild_siblings()
    result=dict(published=PUBLISHED,relationship_names=NAMES,raw_labels=LABELS,n=DISTANCE,d_lin=LINEAL,
                pair_audit=audit,sibling_reconstruction=siblings,
                sources={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [ROOT/'Data/pnas.2300926120.sd01.xlsx',ROOT/'Data/pnas.2300926120.sd03.xlsx']})
    (HERE/'inputs_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(siblings,indent=2));print('Saved inputs_audit.json')

if __name__=='__main__':main()
