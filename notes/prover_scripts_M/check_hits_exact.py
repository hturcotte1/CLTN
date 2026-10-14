"""EXACT follow-up of graphs that are domination-free, legal-stable somewhere (curve test) and have alpha > -1/2:
 (1) adjacency, domination-free check;
 (2) exact alpha test: omega(.;2) Hurwitz minors (not Hurwitz with nonzero Delta_{m-2}, s(2) != 0 => alpha > -1/2);
 (3) curve test cells, and an exact legal stable witness via ctln.analyze_support;
 (4) region2d (repo exact 2-D engine): is the full support an SFP anywhere in the legal plane (positivity included)?"""
import sys, json
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln, graphs, region2d
from curve1d import curve_test
from certify_M_violation import from_d6, hurwitz_minors
from certify_graph_level_pure import omega_at
for s in sys.argv[1:]:
    A = from_d6(s); m = len(A); sig = list(range(m))
    print('=' * 70); print('graph', s, ' A[i][j]=1 iff i->j:')
    for r in A: print('   ', r)
    print(' domination-free:', graphs.is_domination_free(A), ' edges:', sum(map(sum, A)))
    om = omega_at(A, F(2)); Hm = hurwitz_minors(om)
    print(' omega(.;2) Hurwitz minor signs:', [(h > 0) - (h < 0) for h in Hm], ' s(2)=omega(0)=', om[-1], ' hurwitz_stable:', ctln.hurwitz_stable(om)[0])
    r = curve_test(A)
    print(' curve cells:', [(c['lo'][:12], (c['hi'] or 'inf')[:12], c['stable']) for c in r['cells']])
    c = next(c for c in r['cells'] if c['stable'])
    t = F(c['sample_t']); beta = 2 + 1 / t + F(1, 10**6)
    eps = 1 / (1 + t); delta = (beta - 1) / (1 + t)
    rec = ctln.analyze_support(A, sig, eps, delta)
    print(f' witness: t={t} beta={beta} eps={eps} delta={delta} legal={ctln.legal(eps, delta)} stable={rec["stable"]} positive={rec["positive"]}')
    cert = region2d.region2d(A, sig, crosscheck=True)
    print(' region2d verdict:', cert['verdict'], ' positivity cells:', cert.get('positivity_cells'), ' reason:', cert.get('reason'))
    if cert.get('stable_cells'):
        print(' STABLE CELLS:', cert['stable_cells'][:5])
    json.dump(cert, open(f'region2d_{s.replace("&","").replace(chr(92),"_").replace("?","q")}.json', 'w'), default=str)
