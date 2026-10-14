import sys, json, time
sys.path.insert(0, '/home/user/CLTN')
import region2d
rows = ["0011110","0011100","1101111","1110111","0011011","0100101","1100000"]
A = [[int(c) for c in r] for r in rows]
t0 = time.time()
r = region2d.region2d(A, list(range(7)))
r['_elapsed_s'] = round(time.time()-t0, 1)
json.dump(r, open(sys.argv[1], 'w'), indent=1, default=str)
print("verdict", r['verdict'], "stable_cells", r['stable_cells'], "elapsed", r['_elapsed_s'])
