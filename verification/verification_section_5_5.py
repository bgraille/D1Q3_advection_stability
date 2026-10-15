r"""
Verification script for the vertex-stability results (section 5.5)

Parameters
----------
- c   Courant-like parameter
- t_v first rescaled relaxation parameter t_v = s_v - 1;
- t_w second rescaled relaxation parameter t_w = s_w - 1;
- C   Cosine-related quantity used in the symbolic model C = cos(xi); 
- S   Sine-related quantity used in the symbolic model S = sin(xi).

for t_v = -1 and t_w = -1,
the scheme is stable for 0 <= c <= 1

for t_v = +1 and t_w = -1,
the scheme is stable for c = 0 and unstable for 0 < c <= 1

for t_v = -1 and t_w = +1,
the scheme is stable for 0 < c < 1 and unstable for c = 0 or c = 1

for t_v = +1 and t_w = +1,
the scheme is unstable for 0 <= c <= 1

"""
import sage.all as sa    # type: ignore[import-untyped]
import sage.misc.banner  # type: ignore[import-untyped]
from utils import *


sage.misc.banner.SAGE_BANNER = ''
sage.misc.banner.banner()

##############################################################################
# Parameters
##############################################################################
# Management of real and imaginary parts by hands
liste_var = ['c', 't_v', 't_w', 'C', 'S', 'reii', 'imii']
print(f"Create the sage variables {liste_var}")
RingCoeff = sa.PolynomialRing(sa.QQ, names=liste_var)
c, t_v, t_w, C, S, reii, imii = RingCoeff.gens()
XX = sa.var('X', domain=sa.RR)  # indeterminate for polynomials
# Define the ideal associated with the algebraic relations
I = RingCoeff.ideal(
    reii - 1, imii**2 + 1,
    S**2 + C**2 - 1
)
# Define the identity matrix
Id = sa.identity_matrix(RingCoeff, 3)
# Create a polynomial ring with an additional indeterminate lam
RingEigen = sa.PolynomialRing(
    sa.QQ,
    names=["c", "t_v", "t_w", "C", "S", "reii", "imii", "beta"]
)
cE, tv, tw, CE, SE, reiiE, imiiE, beta = RingEigen.gens()
# Define the algebraic relations
IE = RingEigen.ideal(
    reiiE - 1,
    imiiE**2 + 1,
    CE**2+SE**2-1,
    beta**2 - 1 + (1-c**2*(1-CE))**2
)


def ideal_reduce(A, I):
    """Reduce the entries of the matrix A modulo the ideal I."""
    return A.apply_map(lambda a: I.reduce(a))


##############################################################################
# __Proof of Proposition 32__
##############################################################################
banner("Proof of Proposition 32")
print("We assume t_v = -1 and t_w = -1")
valtv, valtw = -1, -1
G = amplification_matrix(c, valtv, valtw, C, S, reii, imii)
print("G = ")
print(G)

##############################################################################
# __Proof of Proposition 33__
##############################################################################
banner("Proof of Proposition 33")
print("We assume t_v = 1 and t_w = -1")
valtv, valtw = 1, -1
G = amplification_matrix(c, valtv, valtw, C, S, reii, imii)
P = G.charpoly(XX)
# Reduce each coefficient of P modulo I
RingPoly = P.parent()
P_reduced = RingPoly([I.reduce(a) for a in P.list()])
print("The characteristic polynomial reads")
print(P_reduced.factor())

print("For c = 0 and 0 <= xi <= pi")
G = amplification_matrix(0, valtv, valtw, C, S, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is X^2-1")
print("check: G^2-Id = ", *ideal_reduce(G*G-Id, I))
assert ideal_reduce(G*G-Id, I).is_zero()

print("For 0 <= c <= 1 and xi = 0")
G = amplification_matrix(c, valtv, valtw, 1, 0, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is X^2-1")
print("check: G^2-Id = ", *ideal_reduce(G*G-Id, I))
assert ideal_reduce(G*G-Id, I).is_zero()

print("For 0 < c <= 1 and xi = pi")
G = amplification_matrix(c, valtv, valtw,-1, 0, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is (X-1)(X^2-1)")
print("check: G^2-Id = ", *ideal_reduce(G*G-Id, I), " != 0")
print("check: (G-Id)(G^2-Id) = ", *ideal_reduce((G-Id)*(G*G-Id), I))
assert ideal_reduce((G-Id)*(G*G-Id), I).is_zero()

print("For c = 1 and xi = pi/2")
G = amplification_matrix(1, valtv, valtw, 0, 1, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is (X-1)(X+i)^2")
print("check: (G-Id)(G+iI) = ", *ideal_reduce((G-Id)*(G+imii*Id), I), " != 0")
print("check: (G-Id)(G+iId)^2 = ", *ideal_reduce((G-Id)*(G+imii*Id)*(G+imii*Id), I))
assert ideal_reduce((G-Id)*(G+imii*Id)*(G+imii*Id), I).is_zero()

##############################################################################
# __Proof of Proposition 34__
##############################################################################
banner("Proof of Proposition 34")
print("We assume t_v = -1 and t_w = 1")
valtv, valtw = -1, 1
G = amplification_matrix(c, valtv, valtw, C, S, reii, imii)
P = G.charpoly(XX)
# Reduce each coefficient of P modulo I
RingPoly = P.parent()
P_reduced = RingPoly([I.reduce(a) for a in P.list()])
print("The characteristic polynomial reads")
print(P_reduced.factor())

print("For c = 0 and 0 < xi <= pi")
G = amplification_matrix(0, valtv, valtw, C, S, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is (X-1)^2(X+1)")
print("check: G^2-Id = ", *ideal_reduce(G*G-Id, I), " != 0")
print("check: (G-Id)^2(G+Id) = ", *ideal_reduce((G+Id)*(G-Id)*(G-Id), I))
assert ideal_reduce((G+Id)*(G-Id)*(G-Id), I).is_zero()

print("For c = 1 and xi=0")
G = amplification_matrix(1, valtv, valtw, 1, 0, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is X^2-1")
print("check: G^2-Id = ", *ideal_reduce(G*G-Id, I))
assert ideal_reduce(G*G-Id, I).is_zero()

print("For c = 1 and xi=pi")
G = amplification_matrix(1, valtv, valtw, -1, 0, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is (X+1)^2")
print("check: G+Id = ", *ideal_reduce(G+Id, I), " != 0")
print("check: (G+Id)^2 = ", *ideal_reduce((G+Id)**2, I))
assert ideal_reduce((G+Id)**2, I).is_zero()

print("For 0 < c < 1 and xi = 0")
G = amplification_matrix(c, valtv, valtw, 1, 0, reii, imii)
print("G = ")
print(ideal_reduce(G, I))
print("The minimal polynomial of G is X^2-1")
print("check: G^2-Id = ", *ideal_reduce(G*G-Id, I))
assert ideal_reduce(G*G-Id, I).is_zero()

print("For 0 < c < 1 and 0 <= xi <= pi")
G = amplification_matrix(c, valtv, valtw, C, S, reii, imii)
print("G = ")
print(ideal_reduce(G, I))

# Transfer G to the new polynomial ring
GE = G.change_ring(RingEigen)
print(GE)
# Matrix associated with the eigenvalue lam
IdE = sa.identity_matrix(RingEigen, 3)
alpha = 1-cE*cE*(1-CE)
lamm = alpha-imiiE*beta
lamp = alpha+imiiE*beta
epixi = CE + imiiE*SE
emixi = CE - imiiE*SE

print("-"*80)
print("Kernel of G - lambda_0 Id (lambda_0 = -1)")
GGE0 = ideal_reduce(GE + IdE, IE)
print(GGE0)
# Extract the first two rows
r1 = sa.vector(RingEigen, GGE0.row(0))
r2 = sa.vector(RingEigen, GGE0.row(1))
# Compute a vector orthogonal to both rows
v0 = r1.cross_product(r2)
# Reduce its components modulo I
v0 = sa.vector(RingEigen, [IE.reduce(a) for a in v0])
print("Candidate eigenvector:")
print(v0)
print("Residual: ", *ideal_reduce(GGE0 * v0.column(), IE))
V0 = sa.vector(RingEigen, [-1-emixi, emixi, 1])
residual = ideal_reduce(GGE0 * V0.column(), IE)
print("Check for ", V0, " : ", *residual)
check_zero(*residual)
print("-"*80)
print("Kernel of G - lambda+ Id (lambda_p± = alpha±isqrt(1-alpha^2) with alpha = 1-c^2(1-C))")
GGEp = (GE - lamp * IdE).apply_map(lambda p: IE.reduce(p))
print(GGEp)
r1 = sa.vector(RingEigen, GGEp.row(0))
r2 = sa.vector(RingEigen, GGEp.row(1))
vp = r1.cross_product(r2)
vp = sa.vector(RingEigen, [IE.reduce(p) for p in vp])
print("Candidate eigenvector:")
print(*vp, sep='\n')
residual = (GGEp * vp.column()).apply_map(lambda p: IE.reduce(p))
print("Residual: ", *residual)
check_zero(*residual)
Vp = sa.vector(
    RingEigen,
    [
        2*(1-cE*cE)*(lamp-emixi),
        (1-cE*cE)*(1-lamp)*emixi,
        -1+lamp-cE*cE*(emixi-lamp*epixi)
    ]
)
residual = ideal_reduce(GGEp * Vp.column(), IE)
print("Check for ", Vp, " : ", *residual)
check_zero(*residual)
print("-"*80)
print("Kernel of G - lambda- Id (lambda_p± = alpha±isqrt(1-alpha^2) with alpha = 1-c^2(1-C))")
GGEm = (GE - lamm * IdE).apply_map(lambda p: IE.reduce(p))
print(GGEm)
r1 = sa.vector(RingEigen, GGEm.row(0))
r2 = sa.vector(RingEigen, GGEm.row(1))
vm = r1.cross_product(r2)
vm = sa.vector(RingEigen, [IE.reduce(p) for p in vm])
print("Candidate eigenvector:")
print(*vm, sep='\n')
residual = (GGEm * vm.column()).apply_map(lambda p: IE.reduce(p))
print("Residual: ", *residual)
check_zero(*residual)
Vm = sa.vector(
    RingEigen,
    [
        2*(1-cE*cE)*(lamm-emixi),
        (1-cE*cE)*(1-lamm)*emixi,
        -1+lamm-cE*cE*(emixi-lamm*epixi)
    ]
)
residual = ideal_reduce(GGEm * Vm.column(), IE)
print("Check for ", Vm, " : ", *residual)
check_zero(*residual)

print("-"*80)
print("Compute the determinant of the change of basis")
P = sa.matrix(RingEigen, [V0, Vm, Vp])
detP = P.det()

print("Matrice de passage P :")
print(P)
print("Déterminant de P :", detP.factor())

##############################################################################
# __Proof of Proposition 35__
##############################################################################
banner("Proof of Proposition 35")
print("We assume t_v = 1 and t_w = 1")
valtv, valtw = 1, 1
epixi = C + imii*S
emixi = C - imii*S

G = amplification_matrix(c, valtv, valtw, C, S, reii, imii)
print(G)
P = G.charpoly(XX)
# Reduce each coefficient of P modulo I
RingPoly = P.parent()
P_reduced = RingPoly([I.reduce(a) for a in P.list()])
print("The characteristic polynomial reads")
print(P_reduced.factor())

P3_re = [
    -t_v*t_w,
    t_v*t_w*(1-c**2*(1-C)) - t_v*((1-C)*c**2+C) - t_w*C,
    -(1-c**2*(1-C)) + t_w*((1-C)*c**2+C) + t_v*C,
    1
]
P3_im = [
    0,
    c*t_w*(1+t_v)*S,
    c*(1+t_v)*S,
    0
]

P1dcs_re, P1dcs_im = CohnSchur(*Derive(P3_re, P3_im))
cond = P1dcs_re[1]**2 + P1dcs_im[1]**2 - P1dcs_re[0]**2 - P1dcs_im[0]**2
cond = simp_substitute(cond, S**2, 1-C**2)
cond = cond.subs({t_v: valtv, t_w: valtw})
print("The condition P3 is a von Neumann polynomial reads")
print(cond.factor())
print(">= 0")
cond = cond // (48*(1-C))
print("<=> xi=0 or ")
print(cond.factor())
print(">= 0")

valc = 0
print("-"*80)
print(f"First case: c={valc}")
print("The condition reads ", cond.subs({c: valc}).factor(), " >= 0")
print("The amplification matrix")
G = amplification_matrix(valc, valtv, valtw, C, S, reii, imii)
print(ideal_reduce(G, I))
print("- For xi=0, the minimal polynomial of G is (X-1)(X+1)")
G = amplification_matrix(valc, valtv, valtw, 1, 0, reii, imii)
print("check: (G-Id)(G+Id) = ", *ideal_reduce(G**2-Id, I))
assert ideal_reduce(G**2-Id, I).is_zero()
print("- For 0 < xi < pi, the minimal polynomial of G is (X-1)(X+e^{-ixi})(X+e^{ixi})")
G = amplification_matrix(valc, valtv, valtw, C, S, reii, imii)
print("check: (G-Id)(G+e^{-ixi})(G+e^{ixi}) = ", *ideal_reduce((G-Id)*(G+emixi*Id)*(G+epixi*Id), I))
assert ideal_reduce((G-Id)*(G+emixi*Id)*(G+epixi*Id), I).is_zero()
print("- For xi=pi, the minimal polynomial of G is (X-1)^2")
G = amplification_matrix(valc, valtv, valtw, -1, 0, reii, imii)
print("check: (G-Id) = ", *ideal_reduce(G-Id, I), " != 0")
print("check: (G-Id)^2 = ", *ideal_reduce((G-Id)**2, I))
assert ideal_reduce((G-Id)**2, I).is_zero()

valc = 1
print("-"*80)
print(f"Second case: c={valc}")
print("The condition reads ", cond.subs({c: valc}).factor(), " >= 0")
print("The amplification matrix")
G = amplification_matrix(valc, valtv, valtw, C, S, reii, imii)
print(ideal_reduce(G, I))
print("- For xi=0, the minimal polynomial of G is (X-1)(X+1)")
G = amplification_matrix(valc, valtv, valtw, 1, 0, reii, imii)
print("check: (G-Id)(G+Id) = ", *ideal_reduce(G**2-Id, I))
assert ideal_reduce(G**2-Id, I).is_zero()
print("- For xi=pi/2, the minimal polynomial of G is (X+1)(X+i)^2")
G = amplification_matrix(valc, valtv, valtw, 0, 1, reii, imii)
print("check: (G+Id)(X+iId) = ", *ideal_reduce((G+Id)*(G+imii*Id), I), " != 0")
print("check: (G+Id)(X+iId)^2 = ", *ideal_reduce((G+Id)*(G+imii*Id)**2, I))
assert ideal_reduce((G+Id)*(G+imii*Id)**2, I).is_zero()
print("- For xi=pi, the minimal polynomial of G is (X+1)^2(X-1)")
G = amplification_matrix(valc, valtv, valtw, -1, 0, reii, imii)
print("check: (G+Id)(G-Id) = ", *ideal_reduce((G+Id)*(G-Id), I), " != 0")
print("check: (G+Id)^2(G-Id) = ", *ideal_reduce((G+Id)**2*(G-Id), I))
assert ideal_reduce((G+Id)**2*(G-Id), I).is_zero()
print("- For 0 < xi < pi/2 or pi/2 < xi < pi, the minimal polynomial of G is (X+1)(X-e^{-ixi})(X+e^{ixi})")
G = amplification_matrix(valc, valtv, valtw, C, S, reii, imii)
print("check: (G+Id)(G-e^{-ixi})(G+e^{ixi}) = ", *ideal_reduce((G+Id)*(G-emixi*Id)*(G+epixi*Id), I))
assert ideal_reduce((G+Id)*(G-emixi*Id)*(G+epixi*Id), I).is_zero()

print("-"*80)
print(f"Third case: 0 < c < 1")
print("- For xi=0, the minimal polynomial of G is (X-1)(X+1)")
G = amplification_matrix(c, valtv, valtw, 1, 0, reii, imii)
print(ideal_reduce(G, I))
print("check: (G-Id)(G+Id) = ", *ideal_reduce(G**2-Id, I))
assert ideal_reduce(G**2-Id, I).is_zero()
print("- For xi=pi, the minimal polynomial of G is (X-1)(X^2-2(1-2c^2)X+1)")
print("The condition reads ", cond.subs({C: -1, S: 0}).factor(), " >= 0")
print("The characteristic polynomial is then a von Neumann polynomial")
G = amplification_matrix(c, valtv, valtw, -1, 0, reii, imii)
print(ideal_reduce(G, I))
print("check: (G-Id)(G^2-2(1-2c^2)G+Id) = ", *ideal_reduce((G-Id)*(G*G-2*(1-2*c*c)*G+Id), I))
assert ideal_reduce((G-Id)*(G*G-2*(1-2*c*c)*G+Id), I).is_zero()
dummy = (2*(1-2*c*c))**2 - 4
print("Discriminant of X^2-2(1-2c^2)X+1: ", dummy.factor(), " < 0 for 0 < c < 1")
print("- For cos(xi)=-1+c^2, the condition reads")
print(cond.subs({C: -1+c*c}).factor())
dummy = -c**4*(1-c**2)*(c**8 + 11*c**4*(1-c**2) + 27*(1-c**2)**2)
print("check: ", dummy - cond.subs({C: -1+c*c}))
check_zero(dummy - cond.subs({C: -1+c*c}))

end_of_verification()
