import sys, random, subprocess
rng = random.Random(11)
def from_d6(s):
    n = ord(s[1]) - 63; p = s[2:]; A = [[0] * n for _ in range(n)]; bit = 0
    for i in range(n):
        for j in range(n):
            c = ord(p[bit // 6]) - 63; A[i][j] = (c >> (5 - bit % 6)) & 1; bit += 1
    return A
def to_d6(A):
    n = len(A); bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(n + 63)
    for k in range(0, len(bits), 6):
        v = 0
        for x in bits[k:k+6]: v = 2*v + x
        s += chr(v + 63)
    return s
lines = [l.split()[0] for l in open(sys.argv[1])]
src = rng.sample(lines, 100000)
# also include the U-dropped sample and the known hits
src += [l.split()[0] for l in open(sys.argv[2])][:5000]
orig, perm = [], []
for s in src:
    A = from_d6(s); n = len(A); assert to_d6(A) == s
    pi = list(range(n)); rng.shuffle(pi)
    B = [[A[pi[i]][pi[j]] for j in range(n)] for i in range(n)]
    orig.append(s); perm.append(to_d6(B))
def run(L):
    out = subprocess.run(['./filter7p_copy'], input='\n'.join(L) + '\n', capture_output=True, text=True).stdout.split('\n')
    return {l.split()[0]: l.split()[2] for l in out if l.strip()}
o, p = run(orig), run(perm)
mism = sum(1 for a, b in zip(orig, perm) if o.get(a, 'drop') != p.get(b, 'drop'))
print('graphs', len(orig), 'kept', len(o), 'mismatch', mism)
