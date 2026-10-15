from sympy import *
q,b,z,x,k,a,lam=symbols('q b z x k a lam')
# sextic b->q
S6b=797*b**6-6468*b**5+22544*b**4-44440*b**3+53280*b**2-38016*b+13824
S6q=797*q**6-1686*q**5+2159*q**4-3004*q**3+2499*q**2-2158*q+1521
print('sextic subst ok', expand(S6b.subs(b,q+1)-S6q)==0)
print('sextic irreducible', factor_list(S6q))
r=[s for s in Poly(S6q).intervals() ]
print('S6q real roots', [N(rt,12) for rt in real_roots(S6q)])
print('561/481=',N(Rational(561,481)),'7/6=',N(Rational(7,6)), 'sign change', S6q.subs(q,Rational(561,481)), S6q.subs(q,Rational(7,6)))
Q5b=b**5-3*b**4+3*b**3-3*b**2+2*b-1
Q5q=q**5+2*q**4+q**3-2*q**2-2*q-1
print('quintic subst ok', expand(Q5b.subs(b,q+1)-Q5q)==0)
print('Q5q roots', [N(rt,10) for rt in real_roots(Q5q)], 'signs at 1,13/12', Q5q.subs(q,1), Q5q.subs(q,Rational(13,12)))
print('Q5b signs at 152/73,25/12', Q5b.subs(b,Rational(152,73)), Q5b.subs(b,Rational(25,12)))
P=q**3-q-1
print('plastic', [N(rt,10) for rt in real_roots(P)], P.subs(q,Rational(45,34)), P.subs(q,Rational(4,3)))
print('b^3-3b^2+2b-1 at q+1', expand((b**3-3*b**2+2*b-1).subs(b,q+1)))
# F_k shift
Rk=(6+k)*z**4+(11+k)*z**3+(16+k)*z**2+10*z+3
sh=expand(16*Rk.subs(z,-Rational(1,2)-x))
print('16Rk shift', collect(sh,x))
c=[sh.coeff(x,4),sh.coeff(x,3),sh.coeff(x,2),sh.coeff(x,1),sh.coeff(x,0)]
def hurwitz(c):
    m=len(c)-1; out=[]
    for kk in range(1,m+1):
        M=Matrix(kk,kk,lambda i,j: c[2*(i+1)-(j+1)] if 0<=2*(i+1)-(j+1)<=m else 0)
        out.append(factor(M.det()))
    return out
print('F hurwitz', hurwitz(c))
# E_a cubic: (a+7) z (z^2+z+1)+2 shifted z=-1/2-y
cc=symbols('cc')
cub=expand((cc*z*(z**2+z+1)+2).subs(z,-Rational(1,2)-x))
cub=expand(-cub)  # make leading positive?
print('E cubic shift', collect(cub,x))
cl=[cub.coeff(x,3),cub.coeff(x,2),cub.coeff(x,1),cub.coeff(x,0)]
print('E hurwitz', hurwitz(cl))
# Theorem A g_k
gk=z**5-3*z**2-2*z-1-a*((k+7)*z**4+(k+12)*z**3+(3*k+21)*z**2+(k+13)*z+(k+5))
h=expand(-gk.subs({a:Rational(200,21)}).subs(z,-Rational(10,21)-x))
hc=[h.coeff(x,i) for i in range(5,-1,-1)]
print('A coeffs', [factor(t) for t in hc])
H=hurwitz(hc)
print('21*Delta1', expand(21*H[0]))
# Prop Mfalse (a): alpha of 011/001/100
A=Matrix([[0,1,1],[0,0,1],[1,0,0]])
V=Matrix([[1,0],[0,1],[-1,-1]])
G=(V.T*V).inv()*V.T*A*V
print('Mfalse charpoly A|perp', factor(G.charpoly(z).as_expr()))
