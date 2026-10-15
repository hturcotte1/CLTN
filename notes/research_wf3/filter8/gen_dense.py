import random, sys
random.seed(int(sys.argv[2])); N = int(sys.argv[1]); n = 8
def d6(A):
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    return "&" + chr(63 + n) + "".join(chr(63 + int("".join(map(str, bits[k:k+6])), 2)) for k in range(0, len(bits), 6))
out = sys.stdout
# controls: K8 and K8 minus 1..3 arcs
K = [[1 if i != j else 0 for j in range(n)] for i in range(n)]
out.write(d6(K) + "\n")
for drop in ([(0,1)], [(0,1),(2,3)], [(0,1),(1,2),(2,0)], [(0,1),(1,0)], [(0,1),(0,2)], [(1,0),(2,0)]):
    B = [r[:] for r in K]
    for (a, b) in drop: B[a][b] = 0
    out.write(d6(B) + "\n")
for _ in range(N):
    pb = random.choice((0.5, 0.65, 0.8)); A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            r = random.random()
            if r < pb: A[i][j] = A[j][i] = 1
            elif r < pb + (1-pb)*0.8:
                if random.random() < 0.5: A[i][j] = 1
                else: A[j][i] = 1
    out.write(d6(A) + "\n")
