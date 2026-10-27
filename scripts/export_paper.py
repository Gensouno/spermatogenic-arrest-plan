#!/usr/bin/env python3
"""Pandoc Markdown-to-Word export; no alternative daily prose source."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from zipfile import ZipFile
import xml.etree.ElementTree as ET
from paper_support import ROOT,read_yaml,frontmatter,digest
from check_paper import check

W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
# Pandoc appends literal w:-prefixed styles to reference styles.xml. Preserve
# that declaration rather than serializing it as ns0 and leaving injected w:
# elements unbound in the resulting Word package.
ET.register_namespace('w',W)
for prefix,namespace in {'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'pic':'http://schemas.openxmlformats.org/drawingml/2006/picture'}.items():ET.register_namespace(prefix,namespace)
def q(name):return '{'+W+'}'+name


def pandoc_path():
    explicit=os.environ.get('PANDOC')
    if explicit:
        if not Path(explicit).is_file():raise FileNotFoundError('PANDOC does not name an existing executable')
        return explicit
    found=shutil.which('pandoc')
    if found:return found
    try:
        import pypandoc
        return pypandoc.get_pandoc_path()
    except (ImportError,OSError) as exc:raise FileNotFoundError('Install local envs/requirements-paper.txt; export never downloads tools automatically') from exc


def rewrite_docx(source,target,reference=False,table_style=True):
    """Generic black academic typography and legible bordered tables; no prose edits."""
    with ZipFile(source) as z:members={n:z.read(n) for n in z.namelist()}
    if reference:
        root=ET.fromstring(members['word/styles.xml'])
        for style in root.findall(q('style')):
            if style.get(q('type'))!='paragraph':continue
            props=style.find(q('rPr'))
            if props is None:props=ET.SubElement(style,q('rPr'))
            colour=props.find(q('color'))
            if colour is None:colour=ET.SubElement(props,q('color'))
            colour.attrib.clear();colour.set(q('val'),'000000')
            font=props.find(q('rFonts'))
            if font is None:font=ET.SubElement(props,q('rFonts'))
            for name in ['ascii','hAnsi','eastAsia','cs']:font.set(q(name),'Times New Roman')
            if style.get(q('styleId'))=='Normal':
                size=props.find(q('sz'))
                if size is None:size=ET.SubElement(props,q('sz'))
                size.set(q('val'),'22')
            pp=style.find(q('pPr'))
            if style.get(q('styleId'))=='TableCaption':
                if pp is None:pp=ET.SubElement(style,q('pPr'))
                if pp.find(q('keepNext')) is None:ET.SubElement(pp,q('keepNext'))
            if pp is not None:
                for border in pp.findall(q('pBdr')):pp.remove(border)
        members['word/styles.xml']=ET.tostring(root,encoding='utf-8',xml_declaration=True)
    doc=ET.fromstring(members['word/document.xml'])
    for section in doc.iter(q('sectPr')) if reference else []:
        size=section.find(q('pgSz'))
        if size is None:size=ET.SubElement(section,q('pgSz'))
        size.set(q('w'),'11906');size.set(q('h'),'16838')
    for table in doc.iter(q('tbl')) if table_style else []:
        props=table.find(q('tblPr'))
        if props is None:props=ET.SubElement(table,q('tblPr'))
        borders=props.find(q('tblBorders'))
        if borders is None:borders=ET.SubElement(props,q('tblBorders'))
        for edge in ['top','left','bottom','right','insideH','insideV']:
            b=borders.find(q(edge))
            if b is None:b=ET.SubElement(borders,q(edge))
            b.attrib.update({q('val'):'single',q('sz'):'4',q('color'):'D9D9D9'})
        for cell in table.iter(q('tc')):
            cp=cell.find(q('tcPr'))
            if cp is None:cp=ET.SubElement(cell,q('tcPr'))
            cb=cp.find(q('tcBorders'))
            if cb is not None:
                for edge in cb:
                    edge.attrib.update({q('val'):'single',q('sz'):'4',q('color'):'D9D9D9'})
        first=table.find(q('tr'))
        rows=table.findall(q('tr'))
        # Keep compact scientific tables together; large tables may paginate
        # normally with their repeating headers.
        if len(rows)<=10 and sum(len(''.join(row.itertext())) for row in rows)<1500:
            for row in rows[:-1]:
                for paragraph in row.iter(q('p')):
                    pp=paragraph.find(q('pPr'))
                    if pp is None:pp=ET.SubElement(paragraph,q('pPr'))
                    if pp.find(q('keepNext')) is None:ET.SubElement(pp,q('keepNext'))
        if first is not None:
            for cell in first.findall(q('tc')):
                cp=cell.find(q('tcPr'))
                if cp is None:cp=ET.SubElement(cell,q('tcPr'))
                shade=cp.find(q('shd'))
                if shade is None:shade=ET.SubElement(cp,q('shd'))
                shade.set(q('fill'),'E8EDF1')
    members['word/document.xml']=ET.tostring(doc,encoding='utf-8',xml_declaration=True)
    with ZipFile(target,'w') as z:
        for name,data in members.items():z.writestr(name,data)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'paper/build/manuscript.docx')
    p.add_argument('--reference-doc',type=Path);args=p.parse_args()
    issues=check(ROOT)
    if any(r['level']=='ERROR' for r in issues):raise SystemExit('Paper check failed; fix evidence/structure before Word export')
    meta=read_yaml(ROOT/'paper/metadata.yaml');pandoc=pandoc_path();output=args.output.resolve()
    if output.suffix!='.docx':raise SystemExit('Output must be .docx')
    if output.exists():raise SystemExit('Output exists; choose a new export path to preserve earlier deliverables')
    output.parent.mkdir(parents=True,exist_ok=True)
    sources=[];hashes={}
    for name in meta['chapters']:
        _,body=frontmatter(ROOT/name)
        sources.append(re.sub(r'<!--.*?-->','',body,flags=re.S));hashes[name]=digest(ROOT/name)
    with tempfile.TemporaryDirectory(prefix='paper-export-') as temp:
        temp=Path(temp);combined=temp/'manuscript.md';combined.write_text('\n\n'.join(sources))
        fields={k:meta[k] for k in ['title','subtitle','lang','keywords'] if meta.get(k)}
        if meta.get('authors'):fields['author']=meta['authors']
        fields['reference-section-title']='References';fields['link-citations']=True
        metadata=temp/'metadata.json';metadata.write_text(json.dumps(fields))
        reference=args.reference_doc
        if reference is None and meta.get('reference_doc'):reference=ROOT/meta['reference_doc']
        custom_reference=reference is not None
        if reference is None:
            original=temp/'default.docx'
            original.write_bytes(subprocess.run([pandoc,'--print-default-data-file','reference.docx'],capture_output=True,check=True).stdout)
            reference=temp/'reference.docx';rewrite_docx(original,reference,reference=True)
        if not reference.is_file():raise SystemExit('Requested reference Word template does not exist')
        raw=temp/'raw.docx'
        command=[pandoc,str(combined),'--from=markdown+tex_math_dollars','--to=docx','--standalone','--citeproc',
            '--bibliography',str(ROOT/'paper/references.bib'),'--metadata-file',str(metadata),
            '--resource-path',os.pathsep.join(str(ROOT/n) for n in ['paper/chapters','paper', '.']),
            '--reference-doc',str(reference),'--output',str(raw)]
        result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,check=True)
        if re.search(r'not found|could not|failed',result.stderr,re.I):raise SystemExit('Pandoc reported unresolved content: '+result.stderr)
        final=temp/'final.docx';rewrite_docx(raw,final,table_style=not custom_reference)
        with ZipFile(final) as z:
            needed={'[Content_Types].xml','word/document.xml','word/styles.xml'}
            if not needed<=set(z.namelist()):raise SystemExit('Invalid Word package')
            for name in z.namelist():
                if name.endswith(('.xml','.rels')):ET.fromstring(z.read(name))
            xml=ET.fromstring(z.read('word/document.xml'));text=''.join(xml.itertext())
            if re.search(r'\[@|<!--\s*claim:',text):raise SystemExit('Unresolved citation/claim marker in Word')
            media=sum(n.startswith('word/media/') for n in z.namelist())
            if '![' in combined.read_text() and not media:raise SystemExit('Missing embedded figures')
        shutil.copy2(final,output)
    report={'output_sha256':digest(output),'chapter_sha256':hashes,'bibliography_sha256':digest(ROOT/'paper/references.bib'),
        'pandoc_version':subprocess.run([pandoc,'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0],
        'embedded_media':media,'structure':'passed','rendered_visual_review':'pending',
        'limits':'OOXML/citations/media checked. Render and inspect every page before delivering; Word fields, fonts and venue styling still need review.'}
    output.with_suffix('.export.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Word export created:',output);print('Visual rendering/review pending; not a final submission.')


if __name__=='__main__':main()
