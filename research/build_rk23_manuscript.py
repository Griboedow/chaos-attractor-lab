"""Build an editable, illustrated manuscript from the frozen RK23 evidence."""
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
from verify_rk23_controller import read,sha
from render_twist import colorize,camera

COLORS=['#007d92','#ad570d','#8c45aa']

def load(root,relative,sources):
    path=root/relative;sources.append(dict(path='results/rk23-controller/'+relative,sha256=sha(path)))
    return json.loads(path.read_text()) if path.suffix=='.json' else read(path)

def finish(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=300,facecolor='white');fig.savefig(out/(name+'.pdf'),facecolor='white');plt.close(fig)

def figures(root,out):
    out.mkdir(parents=True,exist_ok=True);sources=[]
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    hs=load(root,'oscillator/render/portrait-histograms.npz',sources)['histograms'];render=load(root,'oscillator/render/analysis.json',sources)['runs']
    fig,axes=plt.subplots(1,2,figsize=(7,2.5))
    for ax,h,row in zip(axes,hs,render):
        im,_=colorize(h,row['log_density_vmax']);ax.imshow(im,interpolation='nearest',origin='upper');ax.set_axis_off();ax.set_title(row['title'].replace(' / ',' | '),pad=5)
    fig.subplots_adjust(left=.005,right=.995,bottom=.01,top=.88,wspace=.035);finish(fig,out,'figure-1-geometry')
    main=load(root,'oscillator/five-million/follow.npz',sources);hyper=load(root,'oscillator/hyper-five-million/follow.npz',sources)
    for number,ids,folder,long,direction in [(2,[381,354,883],'',main,0),(3,[37,526,689],'hyper-',hyper,1)]:
        fig,axes=plt.subplots(1,2,figsize=(7,2.7));ax=axes[0]
        for ident,color in zip(ids,COLORS):
            rates=long['rates'][list(long['ids']).index(ident)];ax.plot(np.arange(1,21)*.25,rates[:,direction],'-o',ms=2,lw=.8,color=color,label=f'O{ident}')
        ax.axhline(0,color='#7e8c95',lw=.6);ax.set(xlabel='Observed steps (millions)',ylabel=('Leading' if direction==0 else 'Second')+' growth / step',title='Five-million-step check');ax.legend(ncol=3,fontsize=6.5,frameon=False,loc='best',handlelength=1.2,columnspacing=.8)
        ax=axes[1];methods=['70 digits','90 digits','SciPy'] if direction==0 else ['5M','70 digits','90 digits','SciPy']
        for offset,(ident,color) in enumerate(zip(ids,COLORS)):
            series=[]
            if direction==1:series.append(long['rates'][list(long['ids']).index(ident),:,direction])
            for digits in [70,90]:series.append(np.array(load(root,f'oscillator/{folder}precision/case-{ident}-{digits}.json',sources)['rates'])[:,direction])
            series.append(load(root,f'oscillator/{folder}upstream/case-{ident}-analysis.npz',sources)['rates'][:,direction])
            means=np.array([v.mean() for v in series]);lower=np.array([v.min() for v in series]);upper=np.array([v.max() for v in series])
            ax.errorbar(np.arange(len(series))+(offset-1)*.15,means,yerr=[means-lower,upper-means],fmt='o',ms=3,capsize=2,lw=.8,color=color)
        ax.axhline(0,color='#7e8c95',lw=.6);ax.set(xticks=np.arange(len(methods)),xticklabels=methods,title='Independent checks',ylabel=('Leading' if direction==0 else 'Second')+' growth / step')
        fig.subplots_adjust(left=.09,right=.99,bottom=.19,top=.86,wspace=.4);finish(fig,out,f'figure-{number}-evidence')
    fig,axes=plt.subplots(2,3,figsize=(7,3.7))
    for row,ident in enumerate([381,883]):
        d=load(root,f'oscillator/approach/case-{ident}.npz',sources);r=next(v for v in render if v['id']==ident)
        coords=(d['snapshots']-np.array(r['center']))@np.array(r['basis'])@camera(r['azimuth'],r['elevation']).T;half=np.max(np.abs(coords))*1.03
        for col,t in enumerate([0,16,1024]):
            xy=coords[list(d['times']).index(t)];ax=axes[row,col];ax.scatter(xy[:,0],xy[:,1],s=.12,c=COLORS[row*2],alpha=.35,rasterized=True);ax.set(xlim=(-half,half),ylim=(-half,half),aspect='equal',title=f'O{ident} / step {t:,}');ax.set_axis_off()
    fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.92,wspace=.02,hspace=.2);finish(fig,out,'figure-4-approach')
    (out/'provenance.json').write_text(json.dumps(dict(sources=sources,scope='Publication layouts from frozen exact histograms, finite-window arrays and volume snapshots; no new simulation or geometric deformation.'),indent=2),encoding='utf-8')

def main(args):
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/rk23-controller';root=repo/'results/rk23-controller';out=args.output.resolve();out.parent.mkdir(parents=True,exist_ok=True)
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
    scout=read(root/'oscillator/scout/screen.npz');long=read(root/'oscillator/five-million/follow.npz');hyper=read(root/'oscillator/hyper-five-million/follow.npz')
    parameters=[['Case','Damping','Relative tolerance']]+[[f'O{i}',f'{scout["params"][i,0]:.13g}',f'{scout["params"][i,1]:.13g}'] for i in [381,354,883,37,526,689]]
    mainrates=[['Case','Mean full spectrum','Minimum leading window']]+[[f'O{i}',' / '.join(f'{v:.6f}' for v in r.mean(0)),f'{r[:,0].min():.6f}'] for i,r in zip(long['ids'],long['rates'])]
    hrates=[['Minimum second rate','5M','70 digits','90 digits','SciPy']]
    for ident in [37,526,689]:
        vals=[hyper['rates'][list(hyper['ids']).index(ident),:,1].min()]
        for digits in [70,90]:vals.append(np.array(json.loads((root/f'oscillator/hyper-precision/case-{ident}-{digits}.json').read_text())['rates'])[:,1].min())
        vals.append(read(root/f'oscillator/hyper-upstream/case-{ident}-analysis.npz')['rates'][:,1].min());hrates.append([f'O{ident}']+[f'{v:.6f}' for v in vals])
    norm=json.loads((root/'oscillator/error-norm-comparison/report.json').read_text())['runs'];normtable=[['Case','Quotient period','Mean quotient spectrum']]+[[f'O{r["id"]}',r['quotient_period'],' / '.join(f'{v:.6f}' for v in r['means'])] for r in norm]
    tables={
        'parameter_table':(parameters,[.16,.42,.42]),
        'protocol_table':([['Stage','Runs','Warm-up / alignment','Measured steps / windows'],['Scout','1,027','20,000','50,000 / 10'],['Every flag','396','100,000','1,000,000 / 20'],['Selected long checks','6','100,000','5,000,000 / 20'],['Independent precision','12','100,000 / 5,000','100,000 / 10'],['Actual SciPy','6','100,000 / 10,000','290,000 / 29'],['Starts and neighbors','189','100,000','300,000 / 10']],[.27,.11,.27,.35]),
        'main_rate_table':(mainrates,[.13,.57,.3]),'hyper_rate_table':(hrates,[.28,.18,.18,.18,.18]),
        'cohort_table':([['Case','24 solution starts','6 initial steps','32 nearby parameters'],['O381','24 / 24','6 / 6','16 / 32'],['O354','24 / 24','6 / 6','32 / 32'],['O883','24 / 24','6 / 6','24 / 32']],[.16,.28,.25,.31]),
        'norm_table':(normtable,[.16,.26,.58])}
    spec=importlib.util.spec_from_file_location('rk23_content',paper/'content.py');content=importlib.util.module_from_spec(spec);spec.loader.exec_module(content)
    figures_by_kind={'geometry':'figure-1-geometry','main_evidence':'figure-2-evidence','hyper_evidence':'figure-3-evidence','approach':'figure-4-approach'};story=[]
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
        if doc.page>1:canvas.drawString(50,A4[1]-27,'TOLERANCE-SCALE GEOMETRY IN AN ADAPTIVE SOLVER')
        canvas.restoreState()
    title='Tolerance-scale geometry and numerical hyperchaos in an adaptive Runge-Kutta solver'
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=47,bottomMargin=51,title=title,author='Chaos Attractor Lab',subject='Standard RK23: numerical candidates, explicit solution bounds, independent checks and credited prior work')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    writer=PdfWriter();writer.clone_document_from_reader(PdfReader(out))
    for page,data,x,y,scale in VECTOR_DRAWS:writer.pages[page].merge_transformed_page(PdfReader(io.BytesIO(data)).pages[0],Transformation().scale(scale).translate(x,y))
    with out.open('wb') as f:writer.write(f)
    reader=PdfReader(out);record=dict(pages=len(reader.pages),planned_pages=len(content.PAGES),page_text_lengths=[len(p.extract_text()) for p in reader.pages],pdf_sha256=sha(out),content_sha256=sha(paper/'content.py'),builder_sha256=sha(Path(__file__)),evidence_commit='d67a116',scope='Manuscript layout from frozen data. Numerical evidence and exact boundedness remain separate; novelty is unresolved.')
    (out.parent/'build.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
