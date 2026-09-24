"""Audit DOCX chart caches and external workbook links; never follow external links."""
import json
import zipfile
import xml.etree.ElementTree as ET
from build_inputs import ROOT,HERE

NS={'c':'http://schemas.openxmlformats.org/drawingml/2006/chart'}

def main():
    result=[]
    for path in sorted((ROOT/'archive/PNAS_refereeing/ClarkRevise').glob('*.docx')):
        with zipfile.ZipFile(path) as z:
            item=dict(source=str(path.relative_to(ROOT)),embedded_files=[n for n in z.namelist() if '/embeddings/' in n],charts=[],external_links=[])
            for name in z.namelist():
                if name.startswith('word/charts/chart') and name.endswith('.xml'):
                    root=ET.fromstring(z.read(name));series=[]
                    for ser in root.findall('.//c:ser',NS):
                        fields={}
                        for child in ser:
                            tag=child.tag.split('}')[-1]
                            if tag not in ('tx','xVal','yVal','cat','val','errBars'):continue
                            fields[tag]=dict(formulas=[x.text for x in child.findall('.//c:f',NS)],
                                indexed_points=[dict(index=x.get('idx'),value=x.find('c:v',NS).text) for x in child.findall('.//c:pt',NS)],
                                literal_texts=[x.text for x in child.findall('c:v',NS)])
                        series.append(fields)
                    item['charts'].append(dict(part=name,series=series))
                if name.startswith('word/charts/_rels/'):
                    item['external_links'] += [dict(part=name,**el.attrib) for el in ET.fromstring(z.read(name)) if el.get('TargetMode')=='External']
            result.append(item)
    (HERE/'chart_cache_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Saved chart_cache_audit.json; external links not accessed')

if __name__=='__main__':main()
