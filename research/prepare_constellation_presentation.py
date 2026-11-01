"""Register the selected numerical candidates in the common portrait catalog."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'provenance/portrait-catalog.json';catalog=json.loads(path.read_text(encoding='utf-8'))
assert not any(row['id']=='raar-constellations' for row in catalog['rows'])
base='results/raar-constellations/runs/raar-constellation-render01/'
rows=[]
for name,title,ids,featured in [
    ('raar-constellations','RAAR with finite spatial point sets',[1504,1725],True),
    ('raar-constellation-columns','Finite-point RAAR: fragmented columns and sheets',[124,896],False),
    ('raar-constellation-fans','Finite-point RAAR: fans and sparse branches',[611,1873],False)]:
    row=dict(id=name,title=title,group='Numerical algorithms',report='docs/raar-constellations.md',source_portrait=f'docs/assets/experimental/constellation-dense/constellation-{ids[0]}.png',description='Numerical candidates in a known projection algorithm. Five-million-step observations and 33-start cohorts retain positive finite-time growth. Finite-alphabet phase retrieval is established prior art; historical uniqueness remains open.',panels=[dict(label=str(i),source=base+f'case-{i}-histogram.npz',key='histogram',mode='histogram') for i in ids],coordinates='Full three-state iterates, centered and rigidly rotated. Equal spatial units within each panel; independent scales. Two million consecutive states per static portrait.')
    if featured:
        row['featured']=True
        row['animation_panels']=[dict(source=base+f'case-{i}-states.npz',metadata_source=base+f'case-{i}.json',key='states',last_count=250000) for i in ids]
    rows.append(row)
index=next(i for i,r in enumerate(catalog['rows']) if r['id']=='raar')+1
catalog['rows'][index:index]=rows
catalog['rows'][index-1]['description']='Case 1051 has a certified attracting positive-entropy set containing the symbolic dynamics and a Whitney fold. Transitivity and equality with the numerical portrait remain unresolved. Case 1483 is a numerical companion.'
path.write_text(json.dumps(catalog,indent=2)+'\n',encoding='utf-8',newline='\n')
gallery=ROOT/'docs/attractor-gallery.md';text=gallery.read_text(encoding='utf-8')
anchor='### RAAR between spatial circles\n';at=text.index(anchor)
blocks=[]
for row in rows:
    blocks.append(f"### {row['title']}\n\n![{row['title']}](assets/portraits/{row['id']}.png)\n\n{row['description']}\n\n{row['coordinates']}\n\n[Report and prior art](raar-constellations.md)\n\n")
text=text[:at]+''.join(blocks)+text[at:]
old='Case 1051 has a certified symbolic invariant set with entropy at least log(2)/49 per iteration. Attraction of that set is unproved. Case 1483 is a numerical companion.'
text=text.replace(old,catalog['rows'][index-1]['description'])
gallery.write_text(text,encoding='utf-8')
readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8');anchor='### Acceleration memory in a scalar iteration\n'
block='''### Circle and finite-point projections

![Finite-point RAAR, cases 1504 and 1725](docs/assets/portraits/raar-constellations.webp)

[Static portrait](docs/assets/portraits/raar-constellations.png)

**Numerical candidates with layered spatial structure.** A circle and a finite set of points produce separated curved sheets in the three-state RAAR iteration. Cases 1504 and 1725 survive five-million-step observations, 33 initial states and 70/100-digit checks. Both use three sites in the observed post-transient motion. Finite-alphabet phase retrieval and projection-map chaos are known; historical uniqueness of these configurations remains unresolved.

[Model, controls and prior art](docs/raar-constellations.md)

'''
assert anchor in text;text=text.replace(anchor,block+anchor)
anchor='| [Spatial-circle RAAR]';at=text.index(anchor);end=text.index('\n',at)+1
text=text[:end]+'| [Circle / finite-point RAAR](docs/raar-constellations.md) | Layered spatial numerical candidates; three active sites in the selected observations | [Complete survey and controls](results/raar-constellations/README.md), [known phase-retrieval connections](docs/raar-constellations.md#prior-art-and-equivalent-terms) |\n| [Cyclic three-circle RAAR](docs/cyclic-raar-search.md) | Positive finite-time growth, but visually thin configurations | [Complete survey](results/cyclic-raar-search/README.md) |\n'+text[end:]
readme.write_text(text,encoding='utf-8')
style=ROOT/'tools/README-portraits.md';text=style.read_text(encoding='utf-8').replace('The seven main-page portraits','The main-page portraits').replace('states of their archived animation inputs. Motor projections','states of their archived animation inputs. Finite-point RAAR uses the last 250,000\nconsecutive states of each two-million-state input, with the static PCA camera. Motor projections')
style.write_text(text,encoding='utf-8')
print('Registered three presentation rows and the selected rotation.')
