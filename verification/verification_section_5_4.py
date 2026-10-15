r"""
Verification script for the edge-stability results (section 5.4)

Parameters
----------
- c   Courant-like parameter
- t_v first rescaled relaxation parameter t_v = s_v - 1;
- t_w second rescaled relaxation parameter t_w = s_w - 1;
- C   Cosine-related quantity used in the symbolic model C = cos(xi); 
- S   Sine-related quantity used in the symbolic model S = sin(xi).

Considering -1 < t_v < 1 and t_w = -1, 
the lattice Boltzmann scheme is stable for every value 0 <= c <= 1.

Considering t_v = -1 and -1 < t_w < 1,
the lattice Boltzmann scheme is stable for every value 0 <= c <= 1.

Considering -1 < t_v < 1 and t_w = 1,
the lattice Boltzmann scheme is stable for c=0 and c=1
but unstable for every other values 0 < c < 1.

Considering t_v = 1 and -1 < t_w < 1,
the lattice Boltzmann scheme is stable for c=0 and c=1
but unstable for every other values 0 < c < 1.
"""
import sage.all as sa    # type: ignore[import-untyped]
import sage.misc.banner  # type: ignore[import-untyped]
from utils import *


sage.misc.banner.SAGE_BANNER = ''
sage.misc.banner.banner()

##############################################################################
# Parameters
##############################################################################
liste_var = ['c', 't_v', 't_w', 'C', 'S']
print(f"Create the sage variables {liste_var}")
RingCoeff = sa.PolynomialRing(sa.QQ, names=liste_var)
c, t_v, t_w, C, S = RingCoeff.gens()

##############################################################################
# Amplification polynomial
##############################################################################
print("Create the amplification polynomial")
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

##############################################################################
# __Proof of Proposition 28__
##############################################################################
banner("Proof of Proposition 28")
print("We assume -1<t_v<1 and t_w=-1")
valtw = -1

print("First, for 0 <= c < 1 and 0 < xi < pi")
print("the scheme is stable if, and only if P1cscs is a simple von Neumann polynomial")
print("The condition reads")
P1cscs_re, P1cscs_im = CohnSchur(*CohnSchur(P3_re, P3_im))
condcscs = P1cscs_re[1]**2 + P1cscs_im[1]**2 - P1cscs_re[0]**2 - P1cscs_im[0]**2
condcscs = simp_substitute(condcscs, S**2, 1-C**2)
print("psi1cscs = ", condcscs.subs({t_w: valtw}), " >= 0.")
print("Indeed, P1cscs reduces to")
for k in range(len(P1cscs_re)):
    P1cscs_re[k] = P1cscs_re[k].subs({t_w: valtw})
    P1cscs_re[k] = simp_substitute(P1cscs_re[k], S**2, 1-C**2)
    if P1cscs_re[k] != 0:
        P1cscs_re[k] = P1cscs_re[k].factor()
    print(f"Re(alphacscs_{{{1,k}}} = ", P1cscs_re[k])
for k in range(len(P1cscs_im)):
    P1cscs_im[k] = P1cscs_im[k].subs({t_w: valtw})
    P1cscs_im[k] = simp_substitute(P1cscs_im[k], S**2, 1-C**2)
    if P1cscs_im[k] != 0:
        P1cscs_im[k] = P1cscs_im[k].factor()
    print(f"Im(alphacscs_{{{1,k}}} = ", P1cscs_im[k])

print("Second, for c = 1 or C^2 = 1")
print("the scheme is stable if, and only if P1cscde is a simple von Neumann polynomial")
print("The condition reads")
P1csde_re, P1csde_im = Derive(*CohnSchur(P3_re, P3_im))
condcsde = P1csde_re[1]**2 + P1csde_im[1]**2 - P1csde_re[0]**2 - P1csde_im[0]**2
condcsde = simp_substitute(condcsde, S**2, 1-C**2)
print("psi1csde = ", condcsde.subs({t_w: valtw}).factor(), " >= 0.")
print("Indeed, P1csde reduces to")
for k in range(len(P1csde_re)):
    P1csde_re[k] = P1csde_re[k].subs({t_w: valtw})
    P1csde_re[k] = simp_substitute(P1csde_re[k], S**2, 1-C**2)
    if P1csde_re[k] != 0:
        P1csde_re[k] = P1csde_re[k].factor()
    print(f"Re(alphacsde_{{{1,k}}} = ", P1csde_re[k])
for k in range(len(P1csde_im)):
    P1csde_im[k] = P1csde_im[k].subs({t_w: valtw})
    P1csde_im[k] = simp_substitute(P1csde_im[k], S**2, 1-C**2)
    if P1csde_im[k] != 0:
        P1csde_im[k] = P1csde_im[k].factor()
    print(f"Im(alphacsde_{{{1,k}}} = ", P1csde_im[k])

##############################################################################
# __Proof of Proposition 29__
##############################################################################
banner("Proof of Proposition 29")
print("We assume t_v=-1 and -1<t_w<1")
valtv = -1

print("First, for (0 <= c < 1 and 0 < xi < pi) or (0 < c < 1 and xi = pi)")
print("the scheme is stable if, and only if P1cscs is a simple von Neumann polynomial")
print("The condition reads")
P1cscs_re, P1cscs_im = CohnSchur(*CohnSchur(P3_re, P3_im))
condcscs = P1cscs_re[1]**2 + P1cscs_im[1]**2 - P1cscs_re[0]**2 - P1cscs_im[0]**2
condcscs = simp_substitute(condcscs, S**2, 1-C**2)
print("psi1cscs = ", condcscs.subs({t_v: valtv}).factor(), " >= 0.")
print("Indeed, P1cscs reduces to")
for k in range(len(P1cscs_re)):
    P1cscs_re[k] = P1cscs_re[k].subs({t_v: valtv})
    P1cscs_re[k] = simp_substitute(P1cscs_re[k], S**2, 1-C**2)
    if P1cscs_re[k] != 0:
        P1cscs_re[k] = P1cscs_re[k].factor()
    print(f"Re(alphacscs_{{{1,k}}} = ", P1cscs_re[k])
for k in range(len(P1cscs_im)):
    P1cscs_im[k] = P1cscs_im[k].subs({t_v: valtv})
    P1cscs_im[k] = simp_substitute(P1cscs_im[k], S**2, 1-C**2)
    if P1cscs_im[k] != 0:
        P1cscs_im[k] = P1cscs_im[k].factor()
    print(f"Im(alphacscs_{{{1,k}}} = ", P1cscs_im[k])

print("Second, for c = 1 or C = 1 or (c=0 and C=-1)")
print("the scheme is stable if, and only if P1cscde is a simple von Neumann polynomial")
print("The condition reads")
P1csde_re, P1csde_im = Derive(*CohnSchur(P3_re, P3_im))
condcsde = P1csde_re[1]**2 + P1csde_im[1]**2 - P1csde_re[0]**2 - P1csde_im[0]**2
condcsde = simp_substitute(condcsde, S**2, 1-C**2)
print("psi1csde = ", condcsde.subs({t_v: valtv}).factor(), " >= 0.")
print("Indeed, P1csde reduces to")
for k in range(len(P1csde_re)):
    P1csde_re[k] = P1csde_re[k].subs({t_v: valtv})
    P1csde_re[k] = simp_substitute(P1csde_re[k], S**2, 1-C**2)
    if P1csde_re[k] != 0:
        P1csde_re[k] = P1csde_re[k].factor()
    print(f"Re(alphacsde_{{{1,k}}} = ", P1csde_re[k])
for k in range(len(P1csde_im)):
    P1csde_im[k] = P1csde_im[k].subs({t_v: valtv})
    P1csde_im[k] = simp_substitute(P1csde_im[k], S**2, 1-C**2)
    if P1csde_im[k] != 0:
        P1csde_im[k] = P1csde_im[k].factor()
    print(f"Im(alphacsde_{{{1,k}}} = ", P1csde_im[k])

##############################################################################
# __Proof of Proposition 30__
##############################################################################
banner("Proof of Proposition 30")
print("We assume -1<t_v<1 and t_w=1")
valtw = 1

print("First, for 0 <= c < 1 and 0 < xi < pi")
print("the scheme is stable if, and only if P1cscs is a simple von Neumann polynomial")
print("The condition reads")
P1cscs_re, P1cscs_im = CohnSchur(*CohnSchur(P3_re, P3_im))
condcscs = P1cscs_re[1]**2 + P1cscs_im[1]**2 - P1cscs_re[0]**2 - P1cscs_im[0]**2
condcscs = simp_substitute(condcscs, S**2, 1-C**2)
print("psi1cscs = ", condcscs.subs({t_w: valtw}).factor(), " >= 0, if and only if c=0")
print("Indeed, P1cscs reduces to")
for k in range(len(P1cscs_re)):
    P1cscs_re[k] = P1cscs_re[k].subs({t_w: valtw})
    P1cscs_re[k] = simp_substitute(P1cscs_re[k], S**2, 1-C**2)
    if P1cscs_re[k] != 0:
        P1cscs_re[k] = P1cscs_re[k].factor()
    print(f"Re(alphacscs_{{{1,k}}} = ", P1cscs_re[k])
for k in range(len(P1cscs_im)):
    P1cscs_im[k] = P1cscs_im[k].subs({t_w: valtw})
    P1cscs_im[k] = simp_substitute(P1cscs_im[k], S**2, 1-C**2)
    if P1cscs_im[k] != 0:
        P1cscs_im[k] = P1cscs_im[k].factor()
    print(f"Im(alphacscs_{{{1,k}}} = ", P1cscs_im[k])

print("Second, for c = 1 or C^2 = 1")
print("the scheme is stable if, and only if P1cscde is a simple von Neumann polynomial")
print("The condition reads")
P1csde_re, P1csde_im = Derive(*CohnSchur(P3_re, P3_im))
condcsde = P1csde_re[1]**2 + P1csde_im[1]**2 - P1csde_re[0]**2 - P1csde_im[0]**2
condcsde = simp_substitute(condcsde, S**2, 1-C**2)
print("psi1csde = ", condcsde.subs({t_w: valtw}).factor(), " >= 0.")
print("Indeed, P1csde reduces to")
for k in range(len(P1csde_re)):
    P1csde_re[k] = P1csde_re[k].subs({t_w: valtw})
    P1csde_re[k] = simp_substitute(P1csde_re[k], S**2, 1-C**2)
    if P1csde_re[k] != 0:
        P1csde_re[k] = P1csde_re[k].factor()
    print(f"Re(alphacsde_{{{1,k}}} = ", P1csde_re[k])
for k in range(len(P1csde_im)):
    P1csde_im[k] = P1csde_im[k].subs({t_w: valtw})
    P1csde_im[k] = simp_substitute(P1csde_im[k], S**2, 1-C**2)
    if P1csde_im[k] != 0:
        P1csde_im[k] = P1csde_im[k].factor()
    print(f"Im(alphacsde_{{{1,k}}} = ", P1csde_im[k])
print("psi1csde(C=+1) = ", condcsde.subs({t_w: valtw, C: 1, S: 0}).factor(), " >= 0.")
print("psi1csde(C=-1) = ", condcsde.subs({t_w: valtw, C: -1, S: 0}).factor(), " >= 0.")
print("psi1csde(c=+1) = ", condcsde.subs({t_w: valtw, c: 1}).factor(), " >= 0.")

##############################################################################
# __Proof of Proposition 31__
##############################################################################
banner("Proof of Proposition 31")
print("We assume t_v=1 and -1<t_w<1")
valtv = 1

print("First, for (0 <= c < 1 and 0 < xi < pi) or (0 < c < 1 and xi = pi)")
print("the scheme is stable if, and only if P1cscs is a simple von Neumann polynomial")
print("The condition reads")
P1cscs_re, P1cscs_im = CohnSchur(*CohnSchur(P3_re, P3_im))
condcscs = P1cscs_re[1]**2 + P1cscs_im[1]**2 - P1cscs_re[0]**2 - P1cscs_im[0]**2
condcscs = simp_substitute(condcscs, S**2, 1-C**2)
print("psi1cscs = ", condcscs.subs({t_v: valtv}).factor(), " >= 0")
print("psi1cscs >= 0 iff (c=0 and 0<xi<pi) or (0<c<1 and xi=pi).")
print("Indeed, P1cscs reduces to")
for k in range(len(P1cscs_re)):
    P1cscs_re[k] = P1cscs_re[k].subs({t_v: valtv})
    P1cscs_re[k] = simp_substitute(P1cscs_re[k], S**2, 1-C**2)
    if P1cscs_re[k] != 0:
        P1cscs_re[k] = P1cscs_re[k].factor()
    print(f"Re(alphacscs_{{{1,k}}} = ", P1cscs_re[k])
for k in range(len(P1cscs_im)):
    P1cscs_im[k] = P1cscs_im[k].subs({t_v: valtv})
    P1cscs_im[k] = simp_substitute(P1cscs_im[k], S**2, 1-C**2)
    if P1cscs_im[k] != 0:
        P1cscs_im[k] = P1cscs_im[k].factor()
    print(f"Im(alphacscs_{{{1,k}}} = ", P1cscs_im[k])

print("Second, for c = 1 or C = 1 or (c=0 and C=-1)")
print("the scheme is stable if, and only if P1cscde is a simple von Neumann polynomial")
print("The condition reads")
P1csde_re, P1csde_im = Derive(*CohnSchur(P3_re, P3_im))
condcsde = P1csde_re[1]**2 + P1csde_im[1]**2 - P1csde_re[0]**2 - P1csde_im[0]**2
condcsde = simp_substitute(condcsde, S**2, 1-C**2)
print("psi1csde = ", condcsde.subs({t_v: valtv}).factor(), " >= 0.")
print("Indeed, P1csde reduces to")
for k in range(len(P1csde_re)):
    P1csde_re[k] = P1csde_re[k].subs({t_v: valtv})
    P1csde_re[k] = simp_substitute(P1csde_re[k], S**2, 1-C**2)
    if P1csde_re[k] != 0:
        P1csde_re[k] = P1csde_re[k].factor()
    print(f"Re(alphacsde_{{{1,k}}} = ", P1csde_re[k])
for k in range(len(P1csde_im)):
    P1csde_im[k] = P1csde_im[k].subs({t_v: valtv})
    P1csde_im[k] = simp_substitute(P1csde_im[k], S**2, 1-C**2)
    if P1csde_im[k] != 0:
        P1csde_im[k] = P1csde_im[k].factor()
    print(f"Im(alphacsde_{{{1,k}}} = ", P1csde_im[k])
print("psi1csde(C=+1) = ", condcsde.subs({t_w: valtw, C: 1, S: 0}).factor(), " >= 0.")
print("psi1csde(c=+1) = ", condcsde.subs({t_w: valtw, c: 1}).factor(), " >= 0.")
print("psi1csde(C=-1, c=0) = ", condcsde.subs({t_w: valtw, C: -1, S: 0, c: 0}), " >= 0.")

end_of_verification()
