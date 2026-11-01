"""Build the editable spatial-circle RAAR manuscript from frozen evidence."""
import argparse,importlib.util,io,json,re
from pathlib import Path
import numpy as np
from matplotlib.font_manager import findfont,FontProperties
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle
from pypdf import PdfReader,PdfWriter,Transformation
from build_yogi_manuscript import Vector,equation,VECTOR_DRAWS
from build_raar_figures import figures
from verify_rk23_controller import read,sha


def main(args):
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/raar-circles';root=repo/'results/raar-circles';out=args.output.resolve();out.parent.mkdir(parents=True,exist_ok=True)
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
    scout=read(root/'scout/screen.npz');long=read(root/'five-million/follow.npz');cohort=read(root/'cohort/cohort.npz');follow=read(root/'follow/follow.npz')
    parameters=[['Case','Relaxation beta','Radius B','Tilt theta']]+[[f'R{i}']+[f'{v:.12g}' for v in scout['params'][i,:3]] for i in [1051,1483,820]]
    centers=[['Case','Center x','Center y','Center z']]+[[f'R{i}']+[f'{v:.12g}' for v in scout['params'][i,3:]] for i in [1051,1483,820]]
    rates=[['Case','Mean full spectrum / step','Minimum leading window']]
    for i,r in zip(long['ids'],long['rates']):rates.append([f'R{i}',' / '.join(f'{v:.6f}' for v in r.mean(0)),f'{r[:,0].min():.6f}'])
    series=[long['rates'][list(long['ids']).index(820)]]
    for digits in [70,90]:series.append(np.array(json.loads((root/f'precision/case-820-{digits}.json').read_text())['rates']))
    series.append(read(root/'upstream/case-820.npz')['rates'])
    hyper=[['R820 second rate','5M','70 digits','90 digits','Library'],['Mean']+[f'{r[:,1].mean():.6f}' for r in series],['Minimum window']+[f'{r[:,1].min():.6f}' for r in series]]
    fixed=[['Case','Real polynomial roots','Fixed points','Certified result']]
    for r in json.loads((root/'fixed-proof/report.json').read_text())['cases']:fixed.append([f'R{r["id"]}',r['total_real_roots'],r['number_of_fixed_points'],'Both unstable' if r['number_of_fixed_points'] else 'Fixed-point free'])
    approach=[['Case','Median at step 0','Median at step 1,024','Reference sampling median']]
    for i in [1051,1483,820,75]:
        d=read(root/f'approach/case-{i}.npz');approach.append([f'R{i}',f'{np.median(d["distances"][0]):.6g}',f'{np.median(d["distances"][-1]):.6g}',f'{np.median(d["reference_split_distances"]):.6g}'])
    ct=[['Case','32 starts: leading-positive','32 neighbors: leading-positive','Periodic neighbors']]
    for i in long['ids']:
        start=(cohort['ids']==i)&(cohort['groups']=='start');neighbor=(cohort['ids']==i)&(cohort['groups']=='neighbor')
        ct.append([f'R{i}',int(cohort['flags'][start].sum()),int(cohort['flags'][neighbor].sum()),int(np.sum(cohort['period'][neighbor]>0))])
    slices=[['Slice','Settings','Leading-positive','Two-positive','Periodic']]+[[f'R{r["group"]}',r['count'],r['positive'],r['two_positive'],r['periodic']] for r in json.loads((root/'relaxation/report.json').read_text())['groups']]
    tables={'parameter_table':(parameters,[.13,.29,.29,.29]),'center_table':(centers,[.13,.29,.29,.29]),'rate_table':(rates,[.13,.57,.3]),'hyper_table':(hyper,[.28,.18,.18,.18,.18]),'fixed_table':(fixed,[.14,.28,.20,.38]),'approach_table':(approach,[.13,.25,.29,.33]),'cohort_table':(ct,[.13,.29,.32,.26]),'slice_table':(slices,[.16,.18,.24,.22,.20]),'search_table':([['Stage','Settings','Completed','Leading criterion'],['Scout',len(scout['params']),np.sum(scout['counters'][:,0]==0),scout['flags'].sum()],['Every scout flag',len(follow['ids']),np.sum(follow['counters'][:,0]==0),follow['flags'].sum()]],[.3,.18,.22,.3])}
    spec=importlib.util.spec_from_file_location('raar_content',paper/'content.py');content=importlib.util.module_from_spec(spec);spec.loader.exec_module(content)
    figure_names={'geometry':'figure-1-geometry','evidence':'figure-2-evidence','constraints':'figure-3-constraints','approach':'figure-4-approach','relaxation':'figure-5-relaxation','library':'figure-6-library'};story=[]
    def prose(text):
        parts=re.split(r'(<[^>]+>)',text)
        for i,part in enumerate(parts):
            if part.startswith('<'):continue
            for word,char in [('beta','β'),('theta','θ'),('gamma','γ'),('sigma','σ'),('lambda','λ'),('pi','π')]:part=re.sub(r'(?<![A-Za-z])'+word+r'(?![A-Za-z])',char,part)
            parts[i]=re.sub(r'(P|R|c)_([A-Za-z])',r'\1<sub>\2</sub>',part)
        return ''.join(parts)
    for page,(label,blocks) in enumerate(content.PAGES):
        if page:story.append(PageBreak())
        for kind,text in blocks:
            if kind in styles:story.append(Paragraph(prose(text),styles[kind]))
            elif kind=='eq':story.append(equation(text))
            elif kind in tables:story.extend([table(*tables[kind]),Spacer(1,10)])
            elif kind in figure_names:story.append(Vector((figure_dir/(figure_names[kind]+'.pdf')).read_bytes(),width=width,pad=2))
            else:raise ValueError(kind)
    def footer(canvas,doc):
        canvas.saveState();canvas.setStrokeColor(colors.HexColor('#b7c7ce'));canvas.setLineWidth(.4);canvas.line(50,38,A4[0]-50,38);canvas.setFont('Sans',7.2);canvas.setFillColor(colors.HexColor('#536974'));canvas.drawString(50,25,'CHAOS ATTRACTOR LAB / RESEARCH MANUSCRIPT / 2026-11-01');canvas.drawRightString(A4[0]-50,25,str(doc.page))
        if doc.page>1:canvas.drawString(50,A4[1]-27,'RELAXED REFLECTIONS BETWEEN SPATIAL CIRCLES')
        canvas.restoreState()
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=47,bottomMargin=51,title='Serrated arches and layered loops in relaxed reflections between spatial circles',author='Chaos Attractor Lab',subject='Existing RAAR algorithm: numerical candidates, exact geometric certificates and credited prior work')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    writer=PdfWriter();writer.clone_document_from_reader(PdfReader(out))
    for page,data,x,y,scale in VECTOR_DRAWS:writer.pages[page].merge_transformed_page(PdfReader(io.BytesIO(data)).pages[0],Transformation().scale(scale).translate(x,y))
    with out.open('wb') as f:writer.write(f)
    reader=PdfReader(out);record=dict(pages=len(reader.pages),planned_pages=len(content.PAGES),page_text_lengths=[len(p.extract_text()) for p in reader.pages],pdf_sha256=sha(out),content_sha256=sha(paper/'content.py'),builder_sha256=sha(Path(__file__)),figure_builder_sha256=sha(Path(__file__).with_name('build_raar_figures.py')),evidence_manifest_sha256=sha(root/'manifest.json'),scope='Finite numerical evidence, exact case certificates and historical-priority limits remain separate.')
    (out.parent/'build.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
