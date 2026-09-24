"""Create searchable, page-marked source extracts without modifying originals."""
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / 'sources'

def main():
    OUT.mkdir(exist_ok=True)
    paths = list((ROOT / 'Data').glob('*.pdf'))
    archive = ROOT / 'archive/PNAS_refereeing'
    paths += [p for p in archive.rglob('*') if p.suffix.lower() in ('.pdf', '.docx', '.r', '.txt')
              and 'Papers' not in p.parts and p.name not in ('Neale_2003.pdf', 'Genetic Nurture.pdf')]
    paths += [ROOT / 'references/Goldberger_1978_Models_Methods_IQ_Debate_Revised.pdf']
    manifest = []
    for path in sorted(paths):
        rel = path.relative_to(ROOT)
        dest = OUT / ('__'.join(rel.parts) + '.txt')
        if path.suffix.lower() == '.pdf':
            raw = subprocess.check_output(['pdftotext', '-layout', str(path), '-']).decode('utf-8')
            content = '\n'.join(f'\n=== PDF PAGE {i} ===\n{p}' for i,p in enumerate(raw.split('\f'),1) if p.strip())
        elif path.suffix.lower() == '.docx':
            with zipfile.ZipFile(path) as z:
                tree = ET.fromstring(z.read('word/document.xml'))
                ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
                tags = {'{'+ns['w']+'}t','{http://schemas.openxmlformats.org/officeDocument/2006/math}t'}
                content = '\n'.join(''.join(t.text or '' for t in p.iter() if t.tag in tags)
                                    for p in tree.findall('.//w:p',ns))
        else:
            content = path.read_text(encoding='utf-8',errors='replace')
        dest.write_text(f'SOURCE: {rel}\n\n{content}',encoding='utf-8')
        manifest.append(dict(source=str(rel),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),extract=dest.name))
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f'Extracted {len(manifest)} sources to {OUT}')

if __name__ == '__main__':
    main()
