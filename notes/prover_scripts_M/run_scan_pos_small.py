import sys, numpy as np
sys.path.insert(0, '/home/user/CLTN')
from scan_pos import *
from run_scan_small import labeled, iso
taus = np.logspace(-3, 4, 300)
epss = np.linspace(0.002, 0.498, 249)
for legal in [True, False]:
    for m in [2, 3, 4, 5]:
        As = np.array(list(labeled(m)) if m <= 4 else list(iso(m)))
        Ns = N_from_A(As)
        S = W = ST = 0; ex = []
        for b0 in range(0, len(Ns), 1000):
            vs, vw, st = check_pos(Ns[b0:b0+1000], taus, epss, legal=legal)
            S += vs.sum(); W += vw.sum(); ST += st.sum()
            for b in np.nonzero(vs)[0][:3]:
                if len(ex) < 4: ex.append(As[b0+b].astype(int).tolist())
        print('legal' if legal else 'ANY-tau', m, len(As), 'pos&stable-somewhere', ST, 'Ms+viol', S, 'Mw+viol', W)
        for e in ex: print('    ', e)
