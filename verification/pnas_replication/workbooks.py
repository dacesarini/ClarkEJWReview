"""Streaming XLSX reader; no spreadsheet engine or cached formula evaluation."""
from pathlib import PurePosixPath
import zipfile
import xml.etree.ElementTree as ET

NS = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RID = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'

def rows(path, sheet_name):
    """Yield (Excel row number, column-letter -> value) with missing cells absent."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            strings = [''.join(t.text or '' for t in el.findall('.//s:t', NS))
                       for el in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        book = ET.fromstring(z.read('xl/workbook.xml'))
        links = {el.get('Id'):el.get('Target') for el in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheet = next(el for el in book.find('s:sheets', NS) if el.get('name') == sheet_name)
        target = links[sheet.get(RID)]
        target = target.lstrip('/') if target.startswith('/') else str(PurePosixPath('xl')/target)
        with z.open(target) as stream:
            for _,el in ET.iterparse(stream, events=('end',)):
                if el.tag != '{'+NS['s']+'}row':
                    continue
                values = {}
                for c in el:
                    col = ''.join(x for x in c.get('r','') if x.isalpha())
                    v = c.find('s:v', NS)
                    if c.get('t') == 'inlineStr':
                        values[col] = ''.join(t.text or '' for t in c.findall('.//s:t', NS))
                    elif v is not None:
                        values[col] = strings[int(v.text)] if c.get('t')=='s' else v.text
                yield int(el.get('r')), values
                el.clear()

def records(path, sheet):
    iterator = rows(path, sheet)
    _,header = next(iterator)
    if len(set(header.values())) != len(header):
        raise ValueError('Use column-letter rows for duplicate headers')
    for number,row in iterator:
        yield number,{name:row[col] for col,name in header.items() if col in row}
