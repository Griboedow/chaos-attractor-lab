"""Build the folding companion from frozen evidence and editable prose."""
import argparse,importlib.util,io,json,shutil
from decimal import Decimal
from pathlib import Path
from matplotlib.font_manager import findfont,FontProperties
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle
from pypdf import PdfReader,PdfWriter,Transformation
from build_yogi_manuscript import Vector,equation,VECTOR_DRAWS
from verify_rk23_controller import sha


def main(args):
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-folding';root=repo/'results/raar-folding';out=args.output.resolve();out.parent.mkdir(parents=True,exist_ok=True);figure_dir=out.parent/'manuscript-figures';figure_dir.mkdir(exist_ok=True);VECTOR_DRAWS.clear()
    for folder,names in [('render',['raar-cycle-unfolding','raar-cycle-seed-controls']),('inverse-render',['raar-four-preimages'])]:
        for name in names:
            for ext in ['.pdf','.png']:shutil.copy2(root/folder/(name+ext),figure_dir/(name+ext))
    for ext in ['.pdf','.png']:shutil.copy2(repo/'results/raar-phase-classes/final'/('raar-phase-classes'+ext),figure_dir/('raar-phase-classes'+ext))
    for ext in ['.pdf','.png']:shutil.copy2(repo/'results/raar-fold-certificate/render'/('raar-certified-fold'+ext),figure_dir/('raar-certified-fold'+ext))
    for name,family,weight,style in [('Body','DejaVu Serif','normal','normal'),('BodyBold','DejaVu Serif','bold','normal'),('BodyItalic','DejaVu Serif','normal','italic'),('Sans','DejaVu Sans','normal','normal'),('SansBold','DejaVu Sans','bold','normal')]:pdfmetrics.registerFont(TTFont(name,findfont(FontProperties(family=family,weight=weight,style=style))))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBold');pdfmetrics.registerFontFamily('Sans',normal='Sans',bold='SansBold',italic='Sans',boldItalic='SansBold')
    body=ParagraphStyle('body',fontName='Body',fontSize=10,leading=14,spaceAfter=8,textColor=colors.HexColor('#17252c'),allowWidows=0,allowOrphans=0)
    styles={'p':body,'h1':ParagraphStyle('h1',fontName='SansBold',fontSize=17,leading=21,spaceAfter=15,textColor=colors.HexColor('#173b4a'),keepWithNext=True),'h2':ParagraphStyle('h2',fontName='SansBold',fontSize=10.5,leading=14,spaceBefore=6,spaceAfter=7,textColor=colors.HexColor('#173b4a'),keepWithNext=True),'title':ParagraphStyle('title',fontName='SansBold',fontSize=22,leading=27,spaceAfter=12,textColor=colors.HexColor('#173b4a')),'meta':ParagraphStyle('meta',fontName='Sans',fontSize=8,leading=11,spaceAfter=14,textColor=colors.HexColor('#536974')),'caption':ParagraphStyle('caption',parent=body,fontName='Sans',fontSize=8.2,leading=11.2,spaceBefore=5,spaceAfter=11),'ref':ParagraphStyle('ref',parent=body,fontSize=8.4,leading=11.5,spaceAfter=8),'cell':ParagraphStyle('cell',fontName='Sans',fontSize=8.5,leading=11,textColor=colors.HexColor('#17252c'))}
    width=A4[0]-100
    def table(rows,widths):
        t=Table([[Paragraph(str(v),styles['cell']) for v in row] for row in rows],colWidths=[width*v for v in widths],hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eff1')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#94acb7')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#ced9de'))]));return t
    seeds=json.loads((root/'seeds/report.json').read_text())['runs']
    rows=[['Radius','Seeds / completed','Final median reference distance']]+[[format(r['radius'],'.0e'),format(r['completed'],','),format(r['distance_quantiles'][-1][0],'.8f')] for r in seeds]
    tables={'seed_table':(rows,[.2,.3,.5])};figure_names={'unfolding':'raar-cycle-unfolding','controls':'raar-cycle-seed-controls','preimages':'raar-four-preimages','phase':'raar-phase-classes','fold':'raar-certified-fold'}
    spec=importlib.util.spec_from_file_location('folding_content',paper/'content.py');content=importlib.util.module_from_spec(spec);spec.loader.exec_module(content);story=[]
    for j,(label,blocks) in enumerate(content.PAGES):
        if j:story.append(PageBreak())
        for kind,text in blocks:
            if kind in styles:story.append(Paragraph(text,styles[kind]))
            elif kind=='eq':story.append(equation(text))
            elif kind in tables:story.extend([table(*tables[kind]),Spacer(1,9)])
            elif kind in figure_names:story.append(Vector((figure_dir/(figure_names[kind]+'.pdf')).read_bytes(),width=width,pad=2))
            else:raise ValueError(kind)
    def footer(c,doc):
        c.saveState();c.setStrokeColor(colors.HexColor('#b7c7ce'));c.setLineWidth(.4);c.line(50,38,A4[0]-50,38);c.setFont('Sans',7.2);c.setFillColor(colors.HexColor('#536974'));c.drawString(50,25,'CHAOS ATTRACTOR LAB / FOLDING COMPANION / 2026-11-01');c.drawRightString(A4[0]-50,25,str(doc.page))
        if doc.page>1:c.drawString(50,A4[1]-27,'FROM A TINY SEGMENT TO A SERRATED ARCH')
        c.restoreState()
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=47,bottomMargin=51,title='From a tiny segment to a serrated arch: folding in spatial-circle RAAR',author='Chaos Attractor Lab',subject='Finite tangent-seed evolution and complete certified inverse-image counts; explicit priority limits');doc.build(story,onFirstPage=footer,onLaterPages=footer)
    writer=PdfWriter();writer.clone_document_from_reader(PdfReader(out))
    for page,data,x,y,scale in VECTOR_DRAWS:writer.pages[page].merge_transformed_page(PdfReader(io.BytesIO(data)).pages[0],Transformation().scale(scale).translate(x,y))
    with out.open('wb') as f:writer.write(f)
    pages=PdfReader(out).pages;record=dict(pages=len(pages),planned_pages=len(content.PAGES),page_text_lengths=[len(p.extract_text()) for p in pages],pdf_sha256=sha(out),content_sha256=sha(paper/'content.py'),builder_sha256=sha(Path(__file__)),evidence_manifest_sha256=sha(root/'manifest.json'))
    record['additional_manifest_sha256']={name:sha(repo/'results'/name/'manifest.json') for name in ['raar-phase-classes','raar-priority-extension','raar-fold-certificate']}
    (out.parent/'build.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
