"""
enum_cache.py -- enumerate digraphs on n nodes up to isomorphism, in parallel,
and cache the canonical keys to a file so sweep shards never re-enumerate.

Usage:  python3 enum_cache.py --n 6 --procs 3
Output: certificates/digraphs_n<n>.keys   (one 30-bit integer per line, sorted)

Each line encodes the canonical key from graphs.to_key (off-diagonal entries,
row-major) as  sum(bit_k << k).  Decode with key_from_int().

VALIDATION: the count is checked against OEIS A000273 and the run fails loudly
if it does not match. Never use a cache file whose count is wrong.
Same canonical_form as graphs.py -- this only parallelises the augmentation
step of graphs.enumerate_digraphs, it does not change the algorithm.
"""
import argparse, os, sys, time
from multiprocessing import Pool

from graphs import A000273, canonical_form, enumerate_digraphs, from_key

def key_to_int(key):
    v = 0
    for k, b in enumerate(key):
        if b:
            v |= 1 << k
    return v

def key_from_int(v, n):
    m = n * (n - 1)
    return tuple((v >> k) & 1 for k in range(m))

def load_keys(n, path=None):
    path = path or os.path.join("certificates", f"digraphs_n{n}.keys")
    with open(path) as fh:
        ints = [int(line) for line in fh if line.strip()]
    if len(ints) != A000273[n]:
        raise RuntimeError(f"{path}: {len(ints)} keys, expected A000273[{n}]={A000273[n]}")
    return ints

_N = None
def _augment(parent_int):
    n = _N
    base = from_key(key_from_int(parent_int, n - 1), n - 1)
    slots = 2 * (n - 1)
    out = set()
    for mask in range(1 << slots):
        A = [row[:] + [0] for row in base] + [[0] * n]
        for k in range(n - 1):
            if mask >> (2 * k) & 1:
                A[k][n - 1] = 1
            if mask >> (2 * k + 1) & 1:
                A[n - 1][k] = 1
        out.add(key_to_int(canonical_form(A)[0]))
    return out

def _init(n):
    global _N
    _N = n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--out", default="certificates")
    args = ap.parse_args()
    n = args.n
    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, f"digraphs_n{n}.keys")

    t0 = time.time()
    parents = sorted(key_to_int(k) for k in enumerate_digraphs(n - 1))
    assert len(parents) == A000273[n - 1], (len(parents), A000273[n - 1])
    print(f"{len(parents)} parents on {n-1} nodes ({time.time()-t0:.1f}s); augmenting with {args.procs} procs", flush=True)

    seen = set()
    done = 0
    t0 = time.time()
    with Pool(args.procs, initializer=_init, initargs=(n,)) as pool:
        for s in pool.imap_unordered(_augment, parents, chunksize=8):
            seen |= s
            done += 1
            if done % 500 == 0:
                el = time.time() - t0
                print(f"  {done}/{len(parents)} parents  {len(seen)} graphs  {el:.0f}s  eta {el/done*(len(parents)-done):.0f}s", flush=True)
    print(f"augmentation done in {time.time()-t0:.0f}s: {len(seen)} graphs", flush=True)
    if len(seen) != A000273[n]:
        print(f"*** COUNT MISMATCH: got {len(seen)}, A000273[{n}] = {A000273[n]} -- STOP", flush=True)
        sys.exit(1)
    with open(path, "w") as fh:
        for v in sorted(seen):
            fh.write(f"{v}\n")
    print(f"wrote {path}  (count matches A000273[{n}] = {A000273[n]})")

if __name__ == "__main__":
    main()
