"""Build the symbolic-dynamics article with embedded vector figures."""
import argparse,importlib.util,io,json,shutil
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
from archive_raar_symbolic import sha

def main(args):
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-symbolic-dynamics';root=repo/'results/raar-symbolic-dynamics';out=args.output.resolve();out.parent.mkdir(parents=True,exist_ok=True);figs=out.parent/'manuscript-figures';figs.mkdir(exist_ok=True);VECTOR_DRAWS.clear()
    for name in ['geometry','graph','itineraries']:
        for ext in ['pdf','png']:shutil.copy2(root/'runs/raar-symbolic-render02'/f'raar-symbolic-{name}.{ext}',figs/f'raar-symbolic-{name}.{ext}')
    for name,family,weight,style in [('Body','DejaVu Serif','normal','normal'),('BodyBold','DejaVu Serif','bold','normal'),('BodyItalic','DejaVu Serif','normal','italic'),('Sans','DejaVu Sans','normal','normal'),('SansBold','DejaVu Sans','bold','normal')]:pdfmetrics.registerFont(TTFont(name,findfont(FontProperties(family=family,weight=weight,style=style))))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBold');pdfmetrics.registerFontFamily('Sans',normal='Sans',bold='SansBold',italic='Sans',boldItalic='SansBold')
    body=ParagraphStyle('body',fontName='Body',fontSize=10,leading=14,spaceAfter=8,textColor=colors.HexColor('#17252c'),allowWidows=0,allowOrphans=0)
    styles={'p':body,'h1':ParagraphStyle('h1',fontName='SansBold',fontSize=17,leading=21,spaceAfter=15,textColor=colors.HexColor('#173b4a'),keepWithNext=True),'h2':ParagraphStyle('h2',fontName='SansBold',fontSize=10.5,leading=14,spaceBefore=6,spaceAfter=7,textColor=colors.HexColor('#173b4a'),keepWithNext=True),'title':ParagraphStyle('title',fontName='SansBold',fontSize=22,leading=27,spaceAfter=12,textColor=colors.HexColor('#173b4a')),'meta':ParagraphStyle('meta',fontName='Sans',fontSize=8,leading=11,spaceAfter=14,textColor=colors.HexColor('#536974')),'caption':ParagraphStyle('caption',parent=body,fontName='Sans',fontSize=8.2,leading=11.2,spaceBefore=5,spaceAfter=11),'ref':ParagraphStyle('ref',parent=body,fontSize=8.4,leading=11.5,spaceAfter=8),'cell':ParagraphStyle('cell',fontName='Sans',fontSize=8.5,leading=11,textColor=colors.HexColor('#17252c'))}
    width=A4[0]-100
    def table(rows,widths):
        t=Table([[Paragraph(str(v),styles['cell']) for v in row] for row in rows],colWidths=[width*v for v in widths],hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eff1')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#94acb7')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#ced9de'))]));return t
    load=lambda n:json.loads((root/'runs'/n/'report.json').read_text());primary=load('raar-covering-audit-primary01');second=load('raar-covering-audit-ad01');raw=json.loads((root/'runs/raar-covering-refined01/input.json').read_text())
    edges=[['G edge','Primary cells','Second cells']]+[[f'{i} to {j}',f"{sum(r['cells'] for r in primary['checks'] if r['edge']==e):,}",f"{sum(r['cells'] for r in second['checks'] if r['edge']==e):,}"] for e,(i,j) in enumerate(raw['edges'])]+[['Total',f"{primary['total_accepted_cells']:,}",f"{second['total_accepted_cells']:,}"]]
    params=[['Parameter','Exact binary64 input']]+[[n,repr(v)] for n,v in zip(['beta','r','theta','c_x','c_y','c_z'],raw['params'])]
    claims=[['Established here','Not established here'],['Eight covering relations; positive entropy on a compact invariant set','Attraction of that set or equality with the numerical cloud'],['Every binary route sequence has a realization','Unique coding, uniform hyperbolicity or a transverse homoclinic orbit'],['Located least-period-49 saddle with one unstable multiplier','Stability of all longer numerical shooting examples'],['Exact instance and reproducible proof data','Historical uniqueness or a new dynamical equivalence class']]
    tables={'edge_table':(edges,[.3,.35,.35]),'parameter_table':(params,[.25,.75]),'claim_table':(claims,[.5,.5])}
    spec=importlib.util.spec_from_file_location('symbolic_content',paper/'content.py');content=importlib.util.module_from_spec(spec);spec.loader.exec_module(content);story=[]
    for j,(_,blocks) in enumerate(content.PAGES):
        if j:story.append(PageBreak())
        for kind,value in blocks:
            if kind in styles:story.append(Paragraph(value,styles[kind]))
            elif kind=='eq':story.append(equation(value))
            elif kind in tables:story.extend([table(*tables[kind]),Spacer(1,9)])
            elif kind in ['geometry','graph','itineraries']:story.append(Vector((figs/f'raar-symbolic-{kind}.pdf').read_bytes(),width=width,pad=2))
            else:raise ValueError(kind)
    def footer(c,doc):
        c.saveState();c.setStrokeColor(colors.HexColor('#b7c7ce'));c.setLineWidth(.4);c.line(50,38,A4[0]-50,38);c.setFont('Sans',7.2);c.setFillColor(colors.HexColor('#536974'));c.drawString(50,25,'CHAOS ATTRACTOR LAB / SYMBOLIC DYNAMICS / 2026-11-01');c.drawRightString(A4[0]-50,25,str(doc.page))
        if doc.page>1:c.drawString(50,A4[1]-27,'CERTIFIED SYMBOLIC DYNAMICS IN A SPATIAL-CIRCLE RAAR MAP')
        c.restoreState()
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=47,bottomMargin=51,title='Certified symbolic dynamics in a spatial-circle RAAR map',author='Chaos Attractor Lab',subject='Computer-assisted covering relations and positive topological entropy; explicit attraction and priority limits');doc.build(story,onFirstPage=footer,onLaterPages=footer)
    writer=PdfWriter();writer.clone_document_from_reader(PdfReader(out))
    for page,data,x,y,scale in VECTOR_DRAWS:writer.pages[page].merge_transformed_page(PdfReader(io.BytesIO(data)).pages[0],Transformation().scale(scale).translate(x,y))
    with out.open('wb') as f:writer.write(f)
    pages=PdfReader(out).pages;record=dict(pages=len(pages),planned_pages=len(content.PAGES),page_text_lengths=[len(p.extract_text()) for p in pages],pdf_sha256=sha(out),content_sha256=sha(paper/'content.py'),builder_sha256=sha(Path(__file__)),evidence_manifest_sha256=sha(root/'manifest.json'))
    (out.parent/'build.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print(json.dumps(record))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
