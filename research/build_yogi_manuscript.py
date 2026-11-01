"""Build the editable research manuscript and publication-scale figures from frozen evidence."""
import argparse,hashlib,importlib.util,io,json,re,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import mathtext
from matplotlib.font_manager import findfont,FontProperties
from PIL import Image as PILImage
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT,TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Flowable,Image,Table,TableStyle,Preformatted
from pypdf import PdfReader,PdfWriter,Transformation
from render_twist import colorize
from verify_optimizer_followups import arrays,sha
VECTOR_DRAWS=[]
matplotlib.rcParams['mathtext.fontset']='stix'


class Vector(Flowable):
    def __init__(self,data,width=None,pad=5):
        super().__init__()
        self.data=data;box=list(map(float,PdfReader(io.BytesIO(data)).pages[0].mediabox))
        self.native_width=box[2]-box[0];self.native_height=box[3]-box[1];self.target_width=width;self.pad=pad
    def wrap(self,aw,ah):
        self.scale=min(1.,aw/self.native_width) if self.target_width is None else min(self.target_width,aw)/self.native_width
        self.draw_width=self.native_width*self.scale;self.width=aw;self.height=self.native_height*self.scale+2*self.pad
        return self.width,self.height
    def draw(self):
        c=self.canv;x,y=c.absolutePosition((self.width-self.draw_width)/2,self.pad)
        VECTOR_DRAWS.append((c.getPageNumber()-1,self.data,x,y,self.scale))


def equation(text):
    number=re.search(r"\\tag\{([^}]+)\}",text)
    text=re.sub(r"\\tag\{[^}]+\}","",text)
    if number:text+=r"\qquad ("+number.group(1)+")"
    buffer=io.BytesIO();mathtext.math_to_image('$'+text+'$',buffer,prop=FontProperties(family='STIXGeneral',size=11),format='pdf',color='#15252e')
    return Vector(buffer.getvalue(),pad=5)


def figures(repo,out):
    root=repo/'results/yogi-tanh';out.mkdir(parents=True,exist_ok=True);sources=[]
    for f in [root/'render/portrait-histograms.npz',root/'render/analysis.json',root/'five-million/histories.npz',root/'precision/analysis.json',root/'upstream-analysis/windows.npz',root/'upstream-strong-analysis/windows.npz']:
        sources.append(dict(path=f.relative_to(repo).as_posix(),sha256=sha(f)))
    h=arrays(root/'render/portrait-histograms.npz')['histograms'];r=json.loads((root/'render/analysis.json').read_text())['runs']
    fig,axes=plt.subplots(1,2,figsize=(7,2.6),facecolor='white')
    for ax,hist,row in zip(axes,h,r):
        image,_=colorize(hist,row['log_density_vmax']);ax.imshow(image,interpolation='nearest',origin='upper');ax.set_axis_off()
        ax.set_title(f"S{row['id']:04d} / "+('folded ribbons' if row['id']==2850 else 'feathered crescent'),fontsize=10,pad=7,color='#173b4a')
    fig.subplots_adjust(left=.005,right=.995,top=.86,bottom=.01,wspace=.045)
    fig.savefig(out/'figure-1-geometry.png',dpi=300,facecolor='white');fig.savefig(out/'figure-1-geometry.pdf',facecolor='white');plt.close(fig)
    d=arrays(root/'five-million/histories.npz');precision=json.loads((root/'precision/analysis.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(7,4.8),facecolor='white');colorset=['#c36b16','#087c96']
    for row,ident in enumerate([39,7239]):
        w=d['windows'][list(d['ids']).index(ident)]
        ax=axes[row,0]
        for j in [0,1]:ax.plot(np.arange(1,101)*.05,w[:,j],lw=.8,color=colorset[j],label=f'rate {j+1}')
        ax.axhline(0,color='#687885',lw=.6);ax.set_ylim(bottom=-.006);ax.set_xlabel('Observed steps (millions)',fontsize=9);ax.set_ylabel('Log growth / step',fontsize=9);ax.set_title(f'S{ident:04d}: long run',fontsize=10)
        ax.legend(frameon=False,loc='center left',fontsize=8)
        data=[w[:,1]]
        for digits in [70,90]:data.append(np.array(next(t for t in precision['runs'] if t['id']==ident and t['digits']==digits)['windows'])[:,1])
        u=arrays(root/('upstream-analysis' if ident==39 else 'upstream-strong-analysis')/'windows.npz')
        data.append(u['windows'][list(u['ids']).index(ident),:,1]);ax=axes[row,1]
        for j,y in enumerate(data):
            ax.scatter(np.full(len(y),j),y,s=8,alpha=.55,color=colorset[1],linewidths=0)
            ax.plot([j-.2,j+.2],[y.mean(),y.mean()],color=colorset[0],lw=2)
        ax.axhline(0,color='#687885',lw=.6);ax.set_ylim(bottom=-.002);ax.set_xticks(range(4),['5M','70 digits','90 digits','TFA']);ax.set_title('Second rate: every window',fontsize=10)
        for ax in axes[row]:
            ax.tick_params(labelsize=8);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
            for sp in ax.spines.values():sp.set_color('#a9b5bc')
    fig.subplots_adjust(left=.1,right=.98,bottom=.09,top=.94,hspace=.48,wspace=.38)
    fig.savefig(out/'figure-2-evidence.png',dpi=250,facecolor='white');fig.savefig(out/'figure-2-evidence.pdf',facecolor='white');fig.savefig(out/'figure-2-evidence.svg',facecolor='white');plt.close(fig)
    (out/'provenance.json').write_text(json.dumps(dict(sources=sources,scope='No new simulation. Publication layouts read the frozen portrait histograms and exact numerical windows; image origin is upper to retain the recorded camera.'),indent=2)+'\n')


def main(args):
    repo=Path(__file__).resolve().parents[1];paper=repo/'papers/yogi-memory';out=args.output.resolve();out.parent.mkdir(parents=True,exist_ok=True)
    figure_dir=out.parent/(out.stem+'-figures')
    figures(repo,figure_dir)
    for name,family,weight,style in [('Body','DejaVu Serif','normal','normal'),('BodyBold','DejaVu Serif','bold','normal'),('BodyItalic','DejaVu Serif','normal','italic'),('Sans','DejaVu Sans','normal','normal'),('SansBold','DejaVu Sans','bold','normal'),('Mono','DejaVu Sans Mono','normal','normal')]:
        pdfmetrics.registerFont(TTFont(name,findfont(FontProperties(family=family,weight=weight,style=style))))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBold')
    pdfmetrics.registerFontFamily('Sans',normal='Sans',bold='SansBold',italic='Sans',boldItalic='SansBold')
    body=ParagraphStyle('body',fontName='Body',fontSize=10,leading=14,spaceAfter=8,textColor=colors.HexColor('#17252c'),alignment=TA_LEFT,allowWidows=0,allowOrphans=0)
    styles={'p':body,
        'h1':ParagraphStyle('h1',fontName='SansBold',fontSize=17,leading=21,spaceAfter=15,textColor=colors.HexColor('#173b4a'),keepWithNext=True),
        'h2':ParagraphStyle('h2',fontName='SansBold',fontSize=10.5,leading=14,spaceBefore=6,spaceAfter=7,textColor=colors.HexColor('#173b4a'),keepWithNext=True),
        'title':ParagraphStyle('title',fontName='SansBold',fontSize=22,leading=27,spaceAfter=9,textColor=colors.HexColor('#173b4a')),
        'subtitle':ParagraphStyle('subtitle',fontName='Sans',fontSize=13,leading=18,spaceAfter=12,textColor=colors.HexColor('#536974')),
        'meta':ParagraphStyle('meta',fontName='Sans',fontSize=8,leading=11,spaceAfter=14,textColor=colors.HexColor('#536974')),
        'caption':ParagraphStyle('caption',parent=body,fontName='Sans',fontSize=8.2,leading=11.2,spaceBefore=5,spaceAfter=11,alignment=TA_LEFT),
        'ref':ParagraphStyle('ref',parent=body,fontSize=8.9,leading=12.5,spaceAfter=11,alignment=TA_LEFT),
        'cell':ParagraphStyle('cell',fontName='Sans',fontSize=8,leading=11,textColor=colors.HexColor('#17252c')),
        'code':ParagraphStyle('code',fontName='Mono',fontSize=7,leading=10,spaceAfter=12)
    }
    width=A4[0]-100
    def table(rows,widths):
        data=[[Paragraph(str(v),styles['cell']) for v in row] for row in rows]
        tab=Table(data,colWidths=[width*x for x in widths],hAlign='LEFT',repeatRows=1)
        tab.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eff1')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#94acb7')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#ced9de'))]))
        return tab
    tables={
        'parameter_table':([
            ['Case','b','c','e','k'],
            ['S0039','0.9891246411','0.4454789483','6.935501e-5','2446.627400'],
            ['S7239','0.9664813547','0.1086456582','0.001926981471','3341.612641'],
            ['S2850','0.05067918725','0.3904537907','0.001220584372','15.27946360'],
        ],[.12,.22,.22,.22,.22]),
        'protocol_table':([
            ['Stage','Declared runs','Warm-up','Observed / windows'],
            ['Scout','8,256','20,000','50,000 / 5 x 10,000'],
            ['Every scout flag','488','100,000','1,000,000 / 20 x 50,000'],
            ['Long selected checks','8','100,000','5,000,000 / 100 x 50,000'],
            ['Independent precision','6','100,000 + 5,000 alignment','100,000 / 10 x 10,000'],
            ['Fresh-start cohorts','420','100,000','500,000 / 10 x 50,000'],
            ['Actual TensorFlow','10','100,000 + 9,999 transitions','190,000 / 19 x 10,000'],
        ],[.28,.14,.27,.31]),
        'rate_table':([
            ['Minimum second rate','5M binary64','70 digits','90 digits','TensorFlow'],
            ['S0039','0.009550','0.008971','0.008947','0.009087'],
            ['S7239','0.029162','0.027309','0.028068','0.026603'],
        ],[.30,.175,.175,.175,.175]),
        'cohort_table':([
            ['Case / criterion','Autonomous cold','Library-style cold','Nearby parameters','Full-state box'],
            ['S0039 / two','22 / 22','22 / 22','32 / 32','63 / 64'],
            ['S7239 / two','20 / 22','22 / 22','32 / 32','61 / 64'],
            ['S2850 / leading','2 / 22','0 / 22','32 / 32','40 / 64'],
        ],[.25,.1875,.1875,.1875,.1875]),
        'claim_table':([
            ['Claim','Status'],
            ['All legal scalar orbits bounded under Theorem 1 hypotheses','Proved; no historical-priority claim for the technique.'],
            ['Smooth-Yogi symmetric cycle and strict local stability interval','Proved; known Adam family and local equivalence are explicit.'],
            ['Two expanding directions in S0039 and S7239','Supported by finite long-run, independent precision and actual-library checks.'],
            ['New optimizer, first quadratic optimizer chaos or first symmetric cycle','Not claimed; predecessors are credited.'],
            ['Asymptotic chaotic attraction or a new class inequivalent to all prior attractors','Open; not established by this manuscript.'],
        ],[.49,.51])
    }
    spec=importlib.util.spec_from_file_location('manuscript_content',paper/'content.py');content=importlib.util.module_from_spec(spec);spec.loader.exec_module(content)
    story=[]
    def prose(text):
        # Format mathematical prose without changing XML attributes or URLs.
        pieces=re.split(r'(<[^>]+>)',text)
        for i,piece in enumerate(pieces):
            if piece.startswith('<'):continue
            for word,char in [('theta','θ'),('alpha','α'),('sigma','σ'),('kappa','κ')]:
                piece=re.sub(r'\b'+word+r'\b',char,piece)
            for old,new in [('&lt;=',' ≤ '),('&gt;=',' ≥ '),('&lt;',' &lt; '),('&gt;',' &gt; '),('=',' = ')]:
                piece=piece.replace(old,new)
            pieces[i]=re.sub(r' {2,}',' ',piece)
        return ''.join(pieces)
    for page_index,(label,blocks) in enumerate(content.PAGES):
        if page_index:story.append(PageBreak())
        for kind,text in blocks:
            if kind in styles:story.append(Preformatted(text,styles[kind]) if kind=='code' else Paragraph(prose(text),styles[kind]))
            elif kind=='eq':story.append(equation(text))
            elif kind in tables:story.extend([table(*tables[kind]),Spacer(1,12)])
            elif kind in ['portrait','evidence']:
                stem='figure-1-geometry' if kind=='portrait' else 'figure-2-evidence'
                story.append(Vector((figure_dir/f'{stem}.pdf').read_bytes(),width=width,pad=2))
            else:raise ValueError(kind)
    def footer(canvas,doc):
        canvas.saveState();canvas.setStrokeColor(colors.HexColor('#b7c7ce'));canvas.setLineWidth(.4);canvas.line(50,38,A4[0]-50,38)
        canvas.setFont('Sans',7.2);canvas.setFillColor(colors.HexColor('#536974'));canvas.drawString(50,25,'CHAOS ATTRACTOR LAB  /  RESEARCH MANUSCRIPT  /  2026-11-01')
        canvas.drawRightString(A4[0]-50,25,str(doc.page))
        if doc.page>1:canvas.drawString(50,A4[1]-27,'BOUNDED MEMORY AND NUMERICAL HYPERCHAOS')
        canvas.restoreState()
    doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=47,bottomMargin=51,title='Bounded memory and numerical hyperchaos in scalar quadratic optimization',author='Chaos Attractor Lab',subject='Existing smooth Yogi: scalar boundedness, credited local equivalence and finite numerical evidence')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    writer=PdfWriter();writer.clone_document_from_reader(PdfReader(out))
    for page,data,x,y,scale in VECTOR_DRAWS:
        vector=PdfReader(io.BytesIO(data)).pages[0]
        writer.pages[page].merge_transformed_page(vector,Transformation().scale(scale).translate(x,y))
    with out.open('wb') as stream:writer.write(stream)
    reader=PdfReader(out);texts=[p.extract_text() for p in reader.pages]
    assert len(reader.pages)>=len(content.PAGES)
    record=dict(pages=len(reader.pages),planned_sections=len(content.PAGES),page_text_lengths=[len(t) for t in texts],pdf_sha256=sha(out),content_sha256=sha(paper/'content.py'),builder_sha256=sha(Path(__file__)),source_commit='c612e8b',scope='Publication layout from frozen evidence; mathematical arguments and numerical claims remain separate. Visual page review is recorded separately.')
    (out.parent/'build.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
