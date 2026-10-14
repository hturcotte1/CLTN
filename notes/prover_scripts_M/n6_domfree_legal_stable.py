"""HEURISTIC float pass + EXACT confirmation. n=6: domination-free digraphs that are legal-stable somewhere (pure
stability, full support). Float picks candidates; each is then decided EXACTLY by the Lemma-R curve test, and its
alpha(A|perp) vs -1/2 decided exactly (Hurwitz test of omega(.;2)); positivity scanned exactly at the curve's
stable sample points is NOT claimed -- only reported as float."""
import sys, numpy as np, json
sys.path.insert(0, '/home/user/CLTN')
from enum_cache import key_from_int
from graphs import from_key, is_domination_free
from scan_graph_level import nu_inf
from curve1d import curve_test
from certify_graph_level_pure import omega_at
from certify_M_violation import hurwitz_minors
from fractions import Fraction as F
import ctln
from to_d6 import to_d6
cs = 1 + np.logspace(-3, 3.5, 120)
nuP2 = cs - np.sqrt(cs * cs - cs)
keys = [int(l) for l in open('/home/user/CLTN/certificates/digraphs_n6.keys')]
df = []
for v in keys:
    A = from_key(key_from_int(v, 6), 6)
    if is_domination_free(A): df.append(A)
print('domination-free n=6:', len(df), flush=True)
J = np.ones((6, 6)); cand = []
for b0 in range(0, len(df), 20000):
    As = np.array(df[b0:b0 + 20000], float)
    best = np.full(len(As), -np.inf)
    for k, c in enumerate(cs):
        best = np.maximum(best, np.linalg.eigvals(c * J - As).real.min(axis=1) - nuP2[k])
    cand += [df[b0 + i] for i in np.nonzero(best > -1e-3)[0]]
print('float candidates (margin > -1e-3):', len(cand), flush=True)
rows = []
for A in cand:
    r = curve_test(A)
    if r['verdict'] == 'never': continue
    om = omega_at(A, F(2)); hs = ctln.hurwitz_stable(om)[0]
    rows.append({"d6": to_d6(A), "curve": r['verdict'], "omega2_hurwitz": hs,
                 "alpha_lt_minus_half": hs, "stable_cells": [(c['lo'], c['hi']) for c in r['cells'] if c['stable']]})
json.dump(rows, open('n6_domfree_legal_stable.json', 'w'), indent=1)
from collections import Counter
print('exact: domination-free & legal-stable somewhere:', len(rows), Counter((r['curve'], r['omega2_hurwitz']) for r in rows))
