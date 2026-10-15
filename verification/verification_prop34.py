r"""
Verification script for the vertex-stability results (section 5.5)

Parameters
----------
- c   Courant-like parameter
- t_v first rescaled relaxation parameter t_v = s_v - 1;
- t_w second rescaled relaxation parameter t_w = s_w - 1;
- C   Cosine-related quantity used in the symbolic model C = cos(xi); 
- S   Sine-related quantity used in the symbolic model S = sin(xi).
- C2  C2=cos(xi/2)
- S2  S2=sin(xi/2)
- beta formal square root used to define the roots of the amplification polynomial

C + iS = e^(+i xi/2)^2 = (C2+iS2)^2 = C2^2 - S2^2 + 2i C2 S2
beta = sqrt(1-c^2*S2^2) 

In this script, we check the proof of the satbility
for t_v = -1, t_w = +1, and 0 < c < 1

We diagonalize the amplification matrix

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
liste_var = ['c', 't_v', 't_w', 'C2', 'S2', 'reii', 'imii']
print(f"Create the sage variables {liste_var}")
RingCoeff = sa.PolynomialRing(sa.QQ, names=liste_var)
c, t_v, t_w, C2, S2, reii, imii = RingCoeff.gens()
C = C2*C2-S2*S2     # define cos(xi)=cos(xi/2)^2-sin(xi/2)^2
S = 2*C2*S2         # define sin(xi)=2cos(xi/2)sin(xi/2)
XX = sa.var('X', domain=sa.RR)  # indeterminate for polynomials
# Define the ideal associated with the algebraic relations
I = RingCoeff.ideal(
    reii - 1, imii**2 + 1,
    S2**2 + C2**2 - 1,
)
# Define the identity matrix
Id = sa.identity_matrix(RingCoeff, 3)
# Create a polynomial ring with an additional indeterminate lam
RingEigen = sa.PolynomialRing(
    sa.QQ,
    names=["c", "t_v", "t_w", "C2", "S2", "reii", "imii", "beta"]
)
cE, tv, tw, C2E, S2E, reiiE, imiiE, beta = RingEigen.gens()
CE = C2E*C2E-S2E*S2E
SE = 2*C2E*S2E
# Define the algebraic relations
IE = RingEigen.ideal(
    reiiE - 1,
    imiiE**2 + 1,
    C2E**2+S2E**2-1,
    beta**2 - 1+cE**2*S2E**2       # beta = sqrt(1-c^2*sin(xi/2)^2)
)

def ideal_reduce(A, I):
    """Reduce the entries of the matrix A modulo the ideal I."""
    return A.apply_map(lambda a: I.reduce(a))

print("For 0 < c < 1 and 0 <= xi <= pi")
valtv, valtw = -1, 1
G = amplification_matrix(c, valtv, valtw, C, S, reii, imii)
print("G = ")
print(ideal_reduce(G, I))

# Transfer G to the new polynomial ring
GE = G.change_ring(RingEigen)
print(GE)
# Matrix associated with the eigenvalue lam
IdE = sa.identity_matrix(RingEigen, 3)
alpha = 1-cE*cE*(1-CE)
lamm = 1-2*cE**2*S2E**2 - imiiE*2*c*S2*beta #alpha-imiiE*beta
lamp = 1-2*cE**2*S2E**2 + imiiE*2*c*S2*beta #alpha+imiiE*beta
epixi = CE + imiiE*SE
emixi = CE - imiiE*SE
e2pixi = C2E + imiiE*S2E
e2mixi = C2E - imiiE*S2E

##############################################################################
# First eigenvalue: -1
##############################################################################
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
# print(v0)
print(*[v0k.factor() for v0k in v0], sep='\n')
residual = ideal_reduce(GGE0 * v0.column(), IE)
print("Residual: ", *residual)
check_zero(*residual)
V0 = sa.vector(RingEigen, [-1-emixi, emixi, 1])
residual = ideal_reduce(GGE0 * V0.column(), IE)
print("Check for ", V0, " : ", *residual)
check_zero(*residual)
##############################################################################
# Second eigenvalue: lambda_+
##############################################################################
print("-"*80)
print("Kernel of G - lambda+ Id (lambda_p± = alpha±isqrt(1-alpha^2) with alpha = 1-c^2(1-C))")
GGEp = (GE - lamp * IdE).apply_map(lambda p: IE.reduce(p))
print(GGEp)
r1 = sa.vector(RingEigen, GGEp.row(0))
r2 = sa.vector(RingEigen, GGEp.row(1))
vp = r1.cross_product(r2)
vp = sa.vector(RingEigen, [IE.reduce(p) for p in vp])
print("Candidate eigenvector:")
# print(*vp, sep='\n')
print(*[vpk.factor() for vpk in vp], sep='\n')
residual = (GGEp * vp.column()).apply_map(lambda p: IE.reduce(p))
print("Residual: ", *residual)
check_zero(*residual)
Vp = sa.vector(
    RingEigen,
    [
        2*(1-cE*cE),
        cE*e2mixi*(cE*C2E-beta),
        cE*e2pixi*(cE*C2E+beta)
    ]
)
residual = ideal_reduce(GGEp * Vp.column(), IE)
print("Check for ", Vp, " : ", *residual)
check_zero(*residual)
##############################################################################
# Third eigenvalue: lambda_-
##############################################################################
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
        2*(1-cE*cE),
        cE*e2mixi*(cE*C2E+beta),
        cE*e2pixi*(cE*C2E-beta)
    ]
)
residual = ideal_reduce(GGEm * Vm.column(), IE)
print("Check for ", Vm, " : ", *residual)
check_zero(*residual)
print("-"*80)

##############################################################################
# Determinant of the change of basis
##############################################################################
print("Compute the determinant of the change of basis")
P = sa.matrix(RingEigen, [V0, Vm, Vp])
detP = IE.reduce(P.det())

print("Matrice de passage P :")
print(P)
print("Déterminant de P :", detP.factor())

end_of_verification()

