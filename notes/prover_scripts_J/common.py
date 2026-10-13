import sys, json, itertools
sys.path.insert(0, '/home/user/CLTN')
import os
os.chdir('/home/user/CLTN')
import fastregion, graphs, enum_cache, ctln
from fractions import Fraction as F
import sympy as sp

G1 = [[0,1,1,0,1,1],[1,0,1,0,1,1],[1,1,0,1,0,1],[1,1,1,0,0,0],[0,1,1,0,0,0],[0,0,0,1,1,0]]
G2 = [[0,1,1,0,0,0],[0,0,0,1,0,1],[0,0,0,1,1,1],[1,0,0,0,1,1],[1,1,0,1,0,1],[1,1,1,1,1,0]]

_keys = None
def keys():
    global _keys
    if _keys is None:
        _keys = enum_cache.load_keys(6)
    return _keys

def graph_of_idx(idx, n=6):
    return graphs.from_key(enum_cache.key_from_int(keys()[idx], n), n)

def classes():
    d = json.load(open('/home/user/CLTN/certificates/n6_survivor_equivalence_classes.json'))
    return d['classes']

def Nmat(A):
    n = len(A)
    return [[0 if i == j else 1 - A[j][i] for j in range(n)] for i in range(n)]

def Dsig(A):
    bd = fastregion.beta_data(A, list(range(len(A))))
    return (tuple(bd['d']), tuple(tuple(p) for p in bd['phi']), tuple(tuple(p) for p in bd['omega']),
            tuple(sorted(tuple(p) for p in bd['a'])))

def outdeg(A): return [sum(r) for r in A]
def indeg(A): return [sum(A[i][j] for i in range(len(A))) for j in range(len(A))]
def show(A): return ' / '.join(''.join(map(str, r)) for r in A)
