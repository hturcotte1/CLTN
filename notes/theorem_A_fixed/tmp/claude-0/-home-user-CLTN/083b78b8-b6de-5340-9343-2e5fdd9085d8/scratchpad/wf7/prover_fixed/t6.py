import fixlib as fl, sympy as sp
from fractions import Fraction as F
for d6 in ['&F^U|m^QTq?', '&F^vleMoT]?', '&F^vk\\Yo\\U?', '&F^v[\\UpkN?', '&F^v\\eMp{J?', '&F]S`YTbKx?', '&F^TXmOEwd?', '&F^UWxTClb?', '&F^UXk\\Esb?', '&F^UX_TE{h?']:
    A = fl.from_d6(d6); m = len(A)
    D = sp.Matrix(fl.nonedge(A))
    chi = sp.factor((fl.z*sp.eye(m) - D).det())
    rts = [complex(r) for r in sp.Poly(chi, fl.z).nroots()]
    mr = min(r.real for r in rts)
    # exact: all roots Re > -r  <=> Hurwitz of chi(-r - x)
    def hur(rr):
        h = sp.Poly(sp.expand((-1)**m * chi.subs(fl.z, -rr - fl.x)), fl.x)
        return fl.ctln.hurwitz_stable([F(str(sp.Rational(c))) for c in h.all_coeffs()])[0]
    r = sp.Rational(10,21)
    u = (D.T + r*sp.eye(m)).LUsolve(sp.ones(m,1))
    print(d6, chi, 'minRe~%.4f' % mr, 'Hurwitz@-1/2', hur(sp.Rational(1,2)), 'Hurwitz@-10/21', hur(r), 'u>0', all(v > 0 for v in u), 'U', sum(u))
