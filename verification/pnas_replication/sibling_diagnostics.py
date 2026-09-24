"""Investigate sample/orientation sensitivity without choosing rules to match targets."""
from collections import defaultdict
from itertools import combinations
import json
import numpy as np
from build_inputs import ROOT,HERE,COHORTS,binary_stats
from workbooks import records

def main():
    groups=defaultdict(list)
    for _,r in records(ROOT/'Data/pnas.2300926120.sd01.xlsx','Data'):
        if r.get('dfem')!='0' or r.get('ded') not in ('0','1') or 'byr' not in r:continue
        if not r.get('pid_fath') or r['pid_fath']=='0':continue
        if not 1600<=float(r['byr'])<1920:continue
        groups[r['pid_fath']].append(r)
    counts=defaultdict(lambda:np.zeros(4,dtype=int))
    for group in groups.values():
        for a,b in combinations(group,2):
            a,b=sorted([a,b],key=lambda x:(float(x['byr']),float(x['pid'])))
            x,y=int(a['ded']),int(b['ded']);cell=2*x+y
            full=a.get('pid_moth') and a.get('pid_moth')!='0' and a.get('pid_moth')==b.get('pid_moth')
            lineage='general' if all(float(r['nid'])>=1000 for r in (a,b)) else 'elite_or_mixed'
            for cohort,(lo,hi) in COHORTS.items():
                ca=lo<=float(a['byr'])<hi;cb=lo<=float(b['byr'])<hi
                for rule,keep in [('either',ca or cb),('both',ca and cb),('younger',cb),('older',ca)]:
                    if not keep:continue
                    for parent in ['shared_father']+(['same_known_parents'] if full else []):
                        for lin in ['all',lineage]:
                            counts[(cohort,rule,parent,lin)][cell]+=1
    result=[dict(cohort=c,cohort_rule=r,parent_rule=p,lineage=l,orientation='older_birth_first',**binary_stats(v)) for (c,r,p,l),v in sorted(counts.items())]
    (HERE/'sibling_diagnostics.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    for r in result:
        if r['lineage']=='general':print(r['cohort'],r['cohort_rule'],r['parent_rule'],r['N'],round(r['rho'],6))

if __name__=='__main__':main()
