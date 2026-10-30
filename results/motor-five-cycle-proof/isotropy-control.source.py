"""Exact finite-group controls for the five-cluster orbit symmetry corollary."""
import argparse,datetime,hashlib,itertools,json,math
from pathlib import Path

IDENTITY=tuple(range(8))
GROUPS=((0,3),(1,2),(4,),(5,6),(7,))
RELATIVE=(5,0,3,6,7,1,2,4)


def compose(a,b):
    return tuple(a[b[j]] for j in range(8))


def inverse(a):
    return tuple(a.index(j) for j in range(8))


def closure(generators):
    result={IDENTITY};queue=[IDENTITY]
    while queue:
        a=queue.pop()
        for g in generators:
            b=compose(g,a)
            if b not in result:
                result.add(b);queue.append(b)
    return result


def swap(a,b):
    p=list(IDENTITY);p[a],p[b]=p[b],p[a]
    return tuple(p)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output')
    certificate=json.loads(args.certificate.read_text())
    assert certificate['status']=='attracting_five_cluster_cycle_certified'
    assert certificate['krawczyk_passed'] and all(v['passed'] for v in certificate['cluster_distinctions'])
    assert certificate['multiplicities']==[2,2,1,2,1]
    assert certificate['relative_permutation']==[3,0,4,1,2]
    generators=[swap(0,3),swap(1,2),swap(5,6)]
    pointwise=closure(generators)
    normalizer=[];orders={}
    for p in itertools.permutations(range(8)):
        pi=inverse(p)
        if all(compose(compose(p,h),pi) in pointwise for h in generators):
            normalizer.append(p)
            power=IDENTITY
            for n in range(1,97):
                power=compose(p,power)
                if power in pointwise:
                    orders[p]=n;break
            else:
                raise AssertionError('Quotient order not found')
    orbital=closure(generators+[RELATIVE])
    assert len(pointwise)==8 and len(normalizer)==96 and max(orders.values())==6
    assert orders[RELATIVE]==6 and len(orbital)==48
    assert all(g in normalizer for g in orbital)
    cluster_action=[]
    for group in GROUPS:
        image=frozenset(RELATIVE[j] for j in group)
        cluster_action.append(next(i for i,g in enumerate(GROUPS) if frozenset(g)==image))
    assert cluster_action==certificate['relative_permutation']
    remaining=set(itertools.permutations(range(8)));representatives=[]
    while remaining:
        a=min(remaining);coset={compose(h,a) for h in orbital}
        assert len(coset)==48 and coset<=remaining
        representatives.append(a);remaining-=coset
    assert len(representatives)==840 and len(representatives)*len(orbital)==math.factorial(8)
    result=dict(status='passed',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        certificate_path=str(args.certificate),certificate_sha256=hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        pointwise_order=len(pointwise),normalizer_order=len(normalizer),normalizer_quotient_order=len(normalizer)//len(pointwise),
        maximum_cyclic_order_in_normalizer_quotient=max(orders.values()),
        certified_relative_row_permutation=list(RELATIVE),certified_relative_cluster_action=cluster_action,
        relative_quotient_order=orders[RELATIVE],generated_orbital_group_order=len(orbital),labeled_orbit_copies=len(representatives),
        pointwise_group=[list(g) for g in sorted(pointwise)],orbital_group=[list(g) for g in sorted(orbital)],
        left_coset_representatives=[list(g) for g in representatives],
        mathematical_argument='Distinct cluster states give exactly H0=C2^3. Any orbital symmetry normalizes H0 and acts as a phase translation; the kernel is H0. Therefore the finite quotient is cyclic in S3 x S2 and has order at most6. The certified relative closure supplies order6, forcing orbital stabilizer order48. This does not require the displayed period to be least.',
        limitations=['840 labeled symmetry copies of one attracting orbit class; not840 discoveries or a count of every attractor.',
            'The cyclic-phase argument is a mathematical proof stated with these exact finite-group controls, not a proof-assistant derivation.'])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    args.output.with_suffix('.source.py').write_bytes(Path(__file__).read_bytes())
    print('passed',len(normalizer),'normalizer elements',len(orbital),'orbital symmetries',len(representatives),'cosets')


if __name__=='__main__':
    main()
