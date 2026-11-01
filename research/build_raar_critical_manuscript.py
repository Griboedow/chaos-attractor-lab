"""Build the fold-membership companion article from frozen completed evidence."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
import shutil
from matplotlib.font_manager import findfont, FontProperties
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle
from pypdf import PdfReader, PdfWriter, Transformation
from build_yogi_manuscript import Vector, equation, VECTOR_DRAWS
from audit_raar_trapping_batch import file_hash as sha


def main(args):
    repo=Path(__file__).resolve().parents[1];packet=repo/'results/raar-critical-membership';paper=repo/'papers/raar-critical-membership'
    report=json.loads((packet/'runs/raar-critical-membership-audit03/report.json').read_text())
    secondary=json.loads((packet/'runs/raar-critical-membership-audit04/report.json').read_text())
    assert report['status']==secondary['status']=='passed'
    protocol=json.loads((packet/'runs/raar-critical-thin-window01/input.json').read_text())
    out=args.output.resolve();out.parent.mkdir(parents=True,exist_ok=True);figs=out.parent/'manuscript-figures';figs.mkdir(exist_ok=True)
    VECTOR_DRAWS.clear()
    for ext in ['png','pdf']:
        shutil.copy2(packet/'runs/raar-critical-membership-render02'/f'raar-critical-membership.{ext}',figs/f'raar-critical-membership.{ext}')
    for name,family,weight,style in [('Body','DejaVu Serif','normal','normal'),('BodyBold','DejaVu Serif','bold','normal'),('BodyItalic','DejaVu Serif','normal','italic'),('Sans','DejaVu Sans','normal','normal'),('SansBold','DejaVu Sans','bold','normal')]:
        pdfmetrics.registerFont(TTFont(name,findfont(FontProperties(family=family,weight=weight,style=style))))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBold')
    pdfmetrics.registerFontFamily('Sans',normal='Sans',bold='SansBold',italic='Sans',boldItalic='SansBold')
    ink=colors.HexColor('#17252c');accent=colors.HexColor('#173b4a')
    body=ParagraphStyle('body',fontName='Body',fontSize=10,leading=14,spaceAfter=8,textColor=ink,allowWidows=0,allowOrphans=0)
    styles={'p':body,'h1':ParagraphStyle('h1',fontName='SansBold',fontSize=17,leading=21,spaceAfter=15,textColor=accent,keepWithNext=True),
            'h2':ParagraphStyle('h2',fontName='SansBold',fontSize=10.5,leading=14,spaceBefore=6,spaceAfter=7,textColor=accent,keepWithNext=True),
            'title':ParagraphStyle('title',fontName='SansBold',fontSize=22,leading=27,spaceAfter=12,textColor=accent),
            'meta':ParagraphStyle('meta',fontName='Sans',fontSize=8,leading=11,spaceAfter=14,textColor=colors.HexColor('#536974')),
            'caption':ParagraphStyle('caption',parent=body,fontName='Sans',fontSize=8.2,leading=11.2,spaceBefore=5,spaceAfter=11),
            'ref':ParagraphStyle('ref',parent=body,fontSize=8.1,leading=11,spaceAfter=7),
            'cell':ParagraphStyle('cell',fontName='Sans',fontSize=8.5,leading=11,textColor=ink)}
    width=A4[0]-100
    def table(rows,widths):
        t=Table([[Paragraph(str(v),styles['cell']) for v in row] for row in rows],colWidths=[width*v for v in widths],hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eff1')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#94acb7')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#ced9de'))]))
        return t
    params=[['Parameter','Exact binary64 input']]+[[n,repr(v)] for n,v in zip(['beta','r','theta','c_x','c_y','c_z'],protocol['params'])]
    basis=[['B column 1','B column 2','B column 3']]+[[repr(v) for v in row] for row in protocol['bases'][0]]
    stats=[['Complete check','Vector, 70 digits','Scalar AD, 100 digits'],['Self-covering volume leaves','60','41'],['Self-covering face leaves','2','2'],['Opposite-sign face leaves','2','2'],['Regular image cells','452','270'],['Fold-capable image cells','26','13']]
    tables={'parameters':(params,[.25,.75]),'basis':(basis,[1/3,1/3,1/3]),'statistics':(stats,[.50,.25,.25])}
    spec=importlib.util.spec_from_file_location('critical_content',paper/'content.py');content=importlib.util.module_from_spec(spec);spec.loader.exec_module(content)
    story=[]
    for j,(_,blocks) in enumerate(content.PAGES):
        if j:story.append(PageBreak())
        for kind,value in blocks:
            if kind in styles:story.append(Paragraph(value,styles[kind]))
            elif kind=='eq':story.append(equation(value))
            elif kind in tables:story.extend([table(*tables[kind]),Spacer(1,9)])
            elif kind=='figure':story.append(Vector((figs/'raar-critical-membership.pdf').read_bytes(),width=width,pad=2))
            else:raise ValueError(kind)
    def footer(c,doc):
        c.saveState();c.setStrokeColor(colors.HexColor('#b7c7ce'));c.setLineWidth(.4);c.line(50,38,A4[0]-50,38);c.setFont('Sans',7.2);c.setFillColor(colors.HexColor('#536974'));c.drawString(50,25,'CHAOS ATTRACTOR LAB / FOLD MEMBERSHIP / 2026-11-01');c.drawRightString(A4[0]-50,25,str(doc.page))
        if doc.page>1:c.drawString(50,A4[1]-27,'A WHITNEY FOLD IN THE ATTRACTING SPATIAL-CIRCLE RAAR SET')
        c.restoreState()
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=47,bottomMargin=51,title='A Whitney fold in the attracting spatial-circle RAAR set',author='Chaos Attractor Lab',subject='Certified fold membership and ambient local-equivalence obstruction; explicit priority limits')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    writer=PdfWriter();writer.clone_document_from_reader(PdfReader(out))
    for page,data,x,y,scale in VECTOR_DRAWS:
        writer.pages[page].merge_transformed_page(PdfReader(io.BytesIO(data)).pages[0],Transformation().scale(scale).translate(x,y))
    with out.open('wb') as f:writer.write(f)
    pages=PdfReader(out).pages
    record={'pages':len(pages),'planned_pages':len(content.PAGES),'page_text_lengths':[len(p.extract_text()) for p in pages],
            'pdf_sha256':sha(out),'content_sha256':sha(paper/'content.py'),'builder_sha256':sha(Path(__file__)),
            'evidence_manifest_sha256':sha(packet/'manifest.json'),'membership_audit_sha256':sha(packet/'runs/raar-critical-membership-audit03/report.json'),'secondary_audit_sha256':sha(packet/'runs/raar-critical-membership-audit04/report.json'),'theorem_sha256':sha(packet/'theorem.md')}
    (out.parent/'build.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
