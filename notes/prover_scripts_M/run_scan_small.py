import sys, numpy as np, itertools
sys.path.insert(0, '/home/user/CLTN')
from scan_float import *
from enum_cache import key_from_int
from graphs import from_key

taus = np.concatenate([np.logspace(-3, 4, 400)])

def labeled(m):
    pos = [(i, j) for i in range(m) for j in range(m) if i != j]
    for bits in range(1 << len(pos)):
        A = np.zeros((m, m))
        for k, (i, j) in enumerate(pos):
            if bits >> k & 1:
                A[i, j] = 1
        yield A

def iso(m):
    with open(f'/home/user/CLTN/certificates/digraphs_n{m}.keys') as fh:
        for line in fh:
            v = int(line)
            yield np.array(from_key(key_from_int(v, m), m), dtype=float)

if __name__ == "__main__":
  for m in [2, 3, 4, 5]:
      As = list(labeled(m)) if m <= 4 else list(iso(m))
      As = np.array(As)
      Ns = N_from_A(As)
      tot_s = tot_w = 0
      worst = []
      for b0 in range(0, len(Ns), 2000):
          g, gi, vs, vw = check(Ns[b0:b0+2000], taus)
          for b in np.nonzero(vs.any(axis=1))[0]:
              tot_s += 1
              if len(worst) < 5:
                  worst.append((b0 + b, As[b0 + b].astype(int).tolist()))
          tot_w += int(vw.any(axis=1).sum())
      print(m, len(As), 'Ms-viol', tot_s, 'Mw-viol', tot_w)
      for w in worst:
          print('   ', w)
