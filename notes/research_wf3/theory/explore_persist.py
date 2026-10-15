import sympy as sp, itertools, tlib, classext, boundary_tools as bt
A = tlib.from_d6('&FNFL~[udo?')
s3 = sp.sqrt(3)
xv = [sp.Rational(7,2)-s3*sp.I/2, -s3*sp.I, sp.Rational(1,2)+s3*sp.I/2, sp.Rational(1,2)+s3*sp.I/2, -sp.Rational(7,2)-s3*sp.I/2, -1-2*s3*sp.I, 3*s3*sp.I]
Cs = []
for C in itertools.product([0,1], repeat=7):
    s = sp.expand(sum(xv[i] for i in range(7) if C[i]))
    if sp.simplify(s - 1) == 0: Cs.append(C)
print('subsets C with sum x_C = kappa:', Cs)
res = []
for C in Cs:
    T = tuple(1 - c for c in C)
    for S in itertools.product([0,1], repeat=7):
        B = classext.blow(A, S, T, 1)
        r = bt.analyse(B)
        if r['phi3_mult'] == 1 and r['others_strict'] and r.get('R_gt_1/8') and r.get('pos_near_2'):
            res.append((C, S, r))
            print(C, S, r, flush=True)
print(len(res))
