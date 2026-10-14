"""Extra exact checks: (i) 40 random 7-node tournaments: curve_test 'never' via the identically-zero branch, and
e_2 < 0 at 5 random legal points each (ctln.analyze_support stable is False); (ii) all 9,608 5-node digraphs:
every curve-stable graph has B > Nn (the e_2 bound at m = 5), and curve_test == curve_test_fast."""
import sys, json, random, itertools, subprocess
from fractions import Fraction as F
sys.path.insert(0, '/home/user/CLTN')
import ctln
from n7.curve_test import curve_test, curve_test_fast, from_d6, coeff_polys
rng = random.Random(5)
tour = {"n": 0, "never_zero_branch": 0, "points": 0, "stable_true": 0, "zero_coeffs": set()}
for g in range(40):
    A = [[0] * 7 for _ in range(7)]
    for i, j in itertools.combinations(range(7), 2):
        if rng.random() < 0.5: A[i][j] = 1
        else: A[j][i] = 1
    r = curve_test(A); tour["n"] += 1
    if r["verdict"] == "never" and "identically zero" in r.get("reason", ""): tour["never_zero_branch"] += 1; tour["zero_coeffs"].add(r["reason"])
    for k in range(5):
        t = F(rng.randint(1, 10 ** 5), 100); beta = 2 + 1 / t + F(rng.randint(1, 1000), 10 ** rng.randint(3, 9))
        eps = 1 / (t + 1); delta = (beta - 1) * eps; a = ctln.analyze_support(A, list(range(7)), eps, delta)
        tour["points"] += 1; tour["stable_true"] += (a["stable"] is True)
tour["zero_coeffs"] = sorted(tour["zero_coeffs"])
out5 = subprocess.run("/home/user/CLTN/n7/geng -q 5 | /home/user/CLTN/n7/directg -q", shell=True, capture_output=True, text=True).stdout.split()
n5 = {"graphs": len(out5), "stable_somewhere": 0, "never": 0, "escalate": 0, "stable_with_B_le_Nn": 0, "fast_disagree": 0}
for s in out5:
    A = from_d6(s); P = coeff_polys(A); r = curve_test(A, P); rf = curve_test_fast(A)
    n5[r["verdict"]] += 1; n5["fast_disagree"] += (r["verdict"] != rf["verdict"])
    B = sum(1 for i, j in itertools.combinations(range(5), 2) if A[i][j] and A[j][i])
    Nn = sum(1 for i, j in itertools.combinations(range(5), 2) if not A[i][j] and not A[j][i])
    if r["verdict"] == "stable_somewhere" and not B > Nn: n5["stable_with_B_le_Nn"] += 1
res = {"tournaments_n7": tour, "all_n5": n5}
json.dump(res, open('/tmp/claude-0/-home-user-CLTN/083b78b8-b6de-5340-9343-2e5fdd9085d8/scratchpad/wf2/audit_lemmaR/g_extra.json', 'w'), indent=1)
print(json.dumps(res))
