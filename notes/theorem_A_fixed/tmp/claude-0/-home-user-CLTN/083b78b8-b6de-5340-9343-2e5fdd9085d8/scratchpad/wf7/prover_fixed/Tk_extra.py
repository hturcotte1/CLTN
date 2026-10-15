"""
Tk_extra.py -- question (a) for the multi-class families Q2 (classes 1:3, n = 7 + 4k) and Y (classes 4:1, n = 8 + 5k), at fixed
beta beyond R's ceiling (2.8442): exact stable t-sets via Tk_table.stable_a_set (Hurwitz minors linear-in-a family, exact real
root isolation, one exact rational sample per cell).  T_k(beta) := inf{t0 : stable for every legal t > t0}; 'inf' if unstable
for all large t.  Theorem B predicts T_k = legal boundary 1/(beta - 2) for k >= k0 (Q2, beta <= 3: k0 = 6; Y, beta <= 3: k0 = 22;
Y, beta <= 9/2: k0 = 502).  Writes Tk_extra.json.  Decimals are display only.
"""
import json
from Tk_table import run

out = []
out.append(run('Q2', '3', [0, 1, 2, 3, 4, 5, 6, 7, 10, 20]))
out.append(run('Y', '3', [0, 1, 2, 5, 10, 15, 20, 21, 22, 30]))
out.append(run('Y', '9/2', [20, 50, 100, 150, 200, 300, 400, 502]))
json.dump(out, open('Tk_extra.json', 'w'), indent=1)
print('wrote Tk_extra.json')
