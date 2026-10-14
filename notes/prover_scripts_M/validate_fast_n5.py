"""EXACT: curve_test_fast and curve_test give identical verdicts on all 9,608 digraphs with 5 nodes."""
import sys, time
sys.path.insert(0, '/home/user/CLTN')
from curve1d import curve_test_fast, curve_test
from run_scan_small import iso
t0 = time.time(); diff = 0; tally = {}
for A in iso(5):
    A = A.astype(int).tolist()
    v1 = curve_test_fast(A)['verdict']; v2 = curve_test(A)['verdict']
    tally[v1] = tally.get(v1, 0) + 1
    diff += v1 != v2
print('fast tally', tally, 'disagreements', diff, 'sec', round(time.time() - t0))
