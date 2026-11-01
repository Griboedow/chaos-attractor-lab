"""Build an editable, illustrated manuscript from the frozen RK45 evidence."""
import argparse,importlib.util,io,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont,FontProperties
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle
from pypdf import PdfReader,PdfWriter,Transformation
from build_yogi_manuscript import Vector,equation,VECTOR_DRAWS
from verify_rk45_controller import read,sha
from render_twist import colorize,camera

COLORS=['#007d92','#ad570d','#8c45aa']

def load(root,relative,sources):
    path=root/relative;sources.append(dict(path='results/rk45-controller/'+relative,sha256=sha(path)))
    return json.loads(path.read_text()) if path.suffix=='.json' else read(path)

def finish(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=300,facecolor='white');fig.savefig(out/(name+'.pdf'),facecolor='white');plt.close(fig)

from build_rk45_figures import figures


def main(args):
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/rk45-controller';root=repo/'results/rk45-controller';out=args.output.resolve();out.parent.mkdir(parents=True,exist_ok=True)
    figure_dir=out.parent/(out.stem+'-figures');figures(root,figure_dir);VECTOR_DRAWS.clear()
    for name,family,weight,style in [('Body','DejaVu Serif','normal','normal'),('BodyBold','DejaVu Serif','bold','normal'),('BodyItalic','DejaVu Serif','normal','italic'),('Sans','DejaVu Sans','normal','normal'),('SansBold','DejaVu Sans','bold','normal')]:pdfmetrics.registerFont(TTFont(name,findfont(FontProperties(family=family,weight=weight,style=style))))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBold');pdfmetrics.registerFontFamily('Sans',normal='Sans',bold='SansBold',italic='Sans',boldItalic='SansBold')
    body=ParagraphStyle('body',fontName='Body',fontSize=10,leading=14,spaceAfter=8,textColor=colors.HexColor('#17252c'),allowWidows=0,allowOrphans=0)
    styles={'p':body,'h1':ParagraphStyle('h1',fontName='SansBold',fontSize=17,leading=21,spaceAfter=15,textColor=colors.HexColor('#173b4a'),keepWithNext=True),
            'h2':ParagraphStyle('h2',fontName='SansBold',fontSize=10.5,leading=14,spaceBefore=6,spaceAfter=7,textColor=colors.HexColor('#173b4a'),keepWithNext=True),
            'title':ParagraphStyle('title',fontName='SansBold',fontSize=21,leading=26,spaceAfter=9,textColor=colors.HexColor('#173b4a')),
            'subtitle':ParagraphStyle('subtitle',fontName='Sans',fontSize=12,leading=17,spaceAfter=12,textColor=colors.HexColor('#536974')),
            'meta':ParagraphStyle('meta',fontName='Sans',fontSize=8,leading=11,spaceAfter=14,textColor=colors.HexColor('#536974')),
            'caption':ParagraphStyle('caption',parent=body,fontName='Sans',fontSize=8.2,leading=11.2,spaceBefore=5,spaceAfter=11),
            'ref':ParagraphStyle('ref',parent=body,fontSize=8.7,leading=12,spaceAfter=10),
            'cell':ParagraphStyle('cell',fontName='Sans',fontSize=8,leading=11,textColor=colors.HexColor('#17252c'))}
    width=A4[0]-100
    def table(rows,widths):
        tab=Table([[Paragraph(str(v),styles['cell']) for v in row] for row in rows],colWidths=[width*x for x in widths],hAlign='LEFT',repeatRows=1)
        tab.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eff1')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#94acb7')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#ced9de'))]));return tab
    scout=read(root/'oscillator/scout/screen.npz');long=read(root/'oscillator/five-million/follow.npz');cohort=read(root/'oscillator/cohort/cohort.npz');follow=read(root/'oscillator/follow/follow.npz')
    parameters=[['Case','Damping','Relative tolerance']]+[[f'D{i}',f'{scout["params"][i,0]:.13g}',f'{scout["params"][i,1]:.13g}'] for i in [335,448,1890,2006]]
    mainrates=[['Case','Mean full spectrum','Minimum leading window']]
    for ident in [335,448,1890,2006]:
        r=long['rates'][list(long['ids']).index(ident)];mainrates.append([f'D{ident}',' / '.join(f'{v:.6f}' for v in r.mean(0)),f'{r[:,0].min():.6f}'])
    vals=[long['rates'][list(long['ids']).index(335),:,1].min()]
    for digits in [70,90]:vals.append(np.array(json.loads((root/f'oscillator/precision/case-335-{digits}.json').read_text())['rates'])[:,1].min())
    vals.append(read(root/'oscillator/upstream/case-335-analysis.npz')['rates'][:,1].min())
    hrates=[['Minimum second rate','5M','70 digits','90 digits','SciPy'],['D335']+[f'{v:.6f}' for v in vals]]
    norm=json.loads((root/'oscillator/error-norm-comparison/report.json').read_text())['runs'];normtable=[['Case','Quotient period','Mean quotient spectrum']]+[[f'D{r["id"]}',r['quotient_period'],' / '.join(f'{v:.6f}' for v in r['means'])] for r in norm]
    ct=[['Case','24 solution starts','6 initial steps','32 neighbors','Two-positive neighbors']]
    for ident in [335,448,2006,1890]:
        counts=[]
        for group,total in [('start',24),('step',6),('neighbor',32)]:
            mask=(cohort['ids']==ident)&(cohort['groups']==group);counts.append(f'{cohort["flags"][mask].sum()} / {total}')
        mask=(cohort['ids']==ident)&(cohort['groups']=='neighbor');two=np.sum(cohort['flags'][mask]&(cohort['rates'][mask,:,1].min(1)>0));ct.append([f'D{ident}']+counts+[f'{two} / 32'])
    cert=json.loads((root/'source/structure/analysis.json').read_text())['selected_cases'];bt=[['Case','H','K','Squared gain','Bound constant']]
    for r in cert:bt.append([f'D{r["id"]}',r['H'],r['K'],f'{r["squared_gain_decimal"]:.6g}',f'{r["normalized_P_bound_constant"]:.6f}' if r['gain_condition'] else 'Not certified'])
    dr=json.loads((root/'oscillator/damping/report.json').read_text());dt=[['Slice','Settings','Leading positive','Two positive','Periodic']]+[[r['group'],r['count'],r['positive'],r['two_positive'],r['periodic']] for r in dr['groups']]
    tables={
        'parameter_table':(parameters,[.16,.42,.42]),
        'main_rate_table':(mainrates,[.13,.57,.3]),'hyper_rate_table':(hrates,[.28,.18,.18,.18,.18]),
        'cohort_table':(ct,[.12,.23,.19,.20,.26]),'norm_table':(normtable,[.16,.26,.58]),'bound_table':(bt,[.14,.12,.1,.28,.36]),'damping_table':(dt,[.16,.18,.25,.23,.18]),
        'search_table':([['Stage','Settings','Completed','Retained leading criterion'],['Scout',len(scout['params']),np.sum(scout['counters'][:,0]==0),scout['flags'].sum()],['Every scout flag',len(follow['ids']),np.sum(follow['counters'][:,0]==0),follow['flags'].sum()]],[.3,.18,.22,.3])}
    spec=importlib.util.spec_from_file_location('rk45_content',paper/'content.py');content=importlib.util.module_from_spec(spec);spec.loader.exec_module(content)
    figures_by_kind={'geometry':'figure-1-geometry','main_evidence':'figure-2-evidence','approach':'figure-3-approach','damping':'figure-4-damping'};story=[]
    def prose(text):
        parts=re.split(r'(<[^>]+>)',text)
        for i,part in enumerate(parts):
            if part.startswith('<'):continue
            for word,char in [('zeta','ζ'),('rho','ρ'),('tau','τ'),('lambda','λ'),('omega','ω')]:part=re.sub(r'\b'+word+r'\b',char,part)
            parts[i]=part.replace('&lt;=', '≤').replace('&gt;=', '≥').replace('plus or minus', '±')
        return ''.join(parts)
    for page,(label,blocks) in enumerate(content.PAGES):
        if page:story.append(PageBreak())
        for kind,text in blocks:
            if kind in styles:story.append(Paragraph(prose(text),styles[kind]))
            elif kind=='eq':story.append(equation(text))
            elif kind in tables:story.extend([table(*tables[kind]),Spacer(1,10)])
            elif kind in figures_by_kind:story.append(Vector((figure_dir/(figures_by_kind[kind]+'.pdf')).read_bytes(),width=width,pad=2))
            else:raise ValueError(kind)
    def footer(canvas,doc):
        canvas.saveState();canvas.setStrokeColor(colors.HexColor('#b7c7ce'));canvas.setLineWidth(.4);canvas.line(50,38,A4[0]-50,38);canvas.setFont('Sans',7.2);canvas.setFillColor(colors.HexColor('#536974'));canvas.drawString(50,25,'CHAOS ATTRACTOR LAB / RESEARCH MANUSCRIPT / 2026-11-01');canvas.drawRightString(A4[0]-50,25,str(doc.page))
        if doc.page>1:canvas.drawString(50,A4[1]-27,'ROSETTE GEOMETRY IN AN ADAPTIVE DORMAND-PRINCE SOLVER')
        canvas.restoreState()
    title='Twin rosettes and two expanding directions in an adaptive Dormand-Prince solver'
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=47,bottomMargin=51,title=title,author='Chaos Attractor Lab',subject='Standard RK45: numerical candidates, explicit solution bounds, independent checks and credited prior work')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    writer=PdfWriter();writer.clone_document_from_reader(PdfReader(out))
    for page,data,x,y,scale in VECTOR_DRAWS:writer.pages[page].merge_transformed_page(PdfReader(io.BytesIO(data)).pages[0],Transformation().scale(scale).translate(x,y))
    with out.open('wb') as f:writer.write(f)
    reader=PdfReader(out);record=dict(pages=len(reader.pages),planned_pages=len(content.PAGES),page_text_lengths=[len(p.extract_text()) for p in reader.pages],pdf_sha256=sha(out),content_sha256=sha(paper/'content.py'),builder_sha256=sha(Path(__file__)),evidence_manifest_sha256=sha(root/'manifest.json'),scope='Manuscript layout from frozen data. Numerical evidence and exact boundedness remain separate; novelty is unresolved.')
    (out.parent/'build.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
