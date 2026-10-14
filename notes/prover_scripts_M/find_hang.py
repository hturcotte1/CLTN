import sys, signal, time
sys.path.insert(0, '/home/user/CLTN')
from enum_cache import key_from_int
from graphs import from_key, is_domination_free
from curve1d import curve_test
from to_d6 import to_d6
keys = [int(l) for l in open('/home/user/CLTN/certificates/digraphs_n6.keys')]
lo, hi = int(sys.argv[1]), int(sys.argv[2])
class TO(Exception): pass
def h(*a): raise TO()
signal.signal(signal.SIGALRM, h)
for i in range(lo, hi):
    A = from_key(key_from_int(keys[i], 6), 6)
    if not is_domination_free(A): continue
    signal.alarm(20)
    t0 = time.time()
    try:
        curve_test(A)
    except TO:
        print('HANG idx', i, to_d6(A), A, flush=True)
    signal.alarm(0)
print('scan done', lo, hi)
