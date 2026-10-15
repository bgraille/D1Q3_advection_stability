r"""
Verification script for the stability in -1<t_v,t_w<1 (Proposition 27)

Parameters
----------
- c   Courant-like parameter
- t_v first rescaled relaxation parameter t_v = s_v - 1;
- t_w second rescaled relaxation parameter t_w = s_w - 1;
- C   Cosine-related quantity used in the symbolic model C = cos(xi); 
- S   Sine-related quantity used in the symbolic model S = sin(xi).

Assuming -1<t_v<1 and -1<t_w<1, we have
* for c=0 or c=1, the scheme is stable
* for 0<c<1, the scheme is stable if, and only if -1 < t_w <= t_w^*(t_v)

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

print(r"""
In the case -1<t_v,t_w<1, we have proven that
P3 is a simple von Neumann polynomial if, and only if
P3 is a von Neumann polynomial if, and only if
psi1cscs = |\alphacscs_{1,1}|^2-|\alphacscs_{1,0}|^2 >= 0
""")

##############################################################################
# condition psi1cscs = |\alphacscs_{1,1}|^2-|\alphacscs_{1,0}|^2 >= 0
##############################################################################
banner(r"condition psi1cscs = |\alphacscs_{1,1}|^2-|\alphacscs_{1,0}|^2 >= 0")
P1cscs_re, P1cscs_im = CohnSchur(*CohnSchur(P3_re, P3_im))
cond = P1cscs_re[1]**2 + P1cscs_im[1]**2 - P1cscs_re[0]**2 - P1cscs_im[0]**2
cond = simp_substitute(cond, S**2, 1-C**2)
print("psi1cscs = ")
print(cond.factor())

dummy =  c**2 * (1-c**2) * (1-C)**2 * (1 - t_v**2*t_w**2) * (1 + t_w)
theta1cscs = cond // dummy
print("theta1cscs = ")
print(theta1cscs)

print("We write")
print("theta1cscs(xi,c,t_v,t_v) = alpha0 (1-C^2) + alpham/2 (1-C) + alphap/2 (1+C)")
alpha0 = - theta1cscs.coefficient({C: 2}).factor()
print("alpha0 = ")
print(alpha0)
dummy = (c**2*t_v**2*(1-t_w**2)**2 - (t_v+t_w)**2*(1-t_v*t_w)**2)*((1+t_w)*(1-t_v**2)*(1-c**2) + 4*t_v*(1-t_w))
print("check: ", alpha0-dummy)
check_zero(alpha0 - dummy)

alphap = theta1cscs.coefficient({C: 1}).factor() + theta1cscs.coefficient({C: 0}).factor() - alpha0
print(alphap.factor())
alpham = - theta1cscs.coefficient({C: 1}).factor() + theta1cscs.coefficient({C: 0}).factor() - alpha0
print(alpham.factor())
dummy = 2*theta1cscs - 2*alpha0*(1-C**2) - alphap*(1+C) - alpham*(1-C)
print("check: ", dummy)
check_zero(dummy)

print("\nFirst point")
print("-----------\n")
print("For c=0, theta1cscs = ")
condposcz = theta1cscs.subs({c: 0})
print(condposcz.factor())
gammabarucscsz = 1 + t_v*t_w
gammabarucscsu = - t_v - t_w
gammaucscsz = 1-3*t_w + t_v**3*t_w*(3-t_w) + t_v*(t_w**2+t_w-4) - t_v**2*(1+t_w-4*t_w**2)
gammaucscsu = - (t_v+t_w)*((1-t_v**2)*(1+t_w) + 4*t_v*(1-t_w))
dummy = condposcz - (1-t_v*t_w)**2*(gammabarucscsz + gammabarucscsu*C)*(gammaucscsz+gammaucscsu*C)
print("check: ", dummy)
check_zero(dummy)

dummy = gammabarucscsz - gammabarucscsu
print("gammabarucscsz - gammabarucscsu = ", dummy.factor(), " > 0")
dummy = gammabarucscsz + gammabarucscsu
print("gammabarucscsz + gammabarucscsu = ", dummy.factor(), " > 0")
dummy = gammaucscsz - gammaucscsu
print("gammaucscsz - gammaucscsu = ", dummy.factor(), " > 0")
dummy = gammaucscsz + gammaucscsu
print("gammaucscsz + gammaucscsu = ")
print(dummy.factor())

print("\nSecond point")
print("------------\n")
print("We assume that -1<t_v<1 and -1<t_w<t_w^*")

print("theta1cscs(C=+1) = ", theta1cscs.subs({C: 1}).factor(), " > 0")
print("theta1cscs(C=-1) = ", theta1cscs.subs({C:-1}).factor(), " > 0")

print(r"Compute Delta the discriminant of the polynomial theta1cscs in the variable C")
print("We have to prove that this discriminant is always non negative "
      "in order to prove that the roots are always real.\n"
      "We will then avoid complex roots that go inside the square "
      "0 <= c <= 1, -1 <= C <= 1.")
cA = theta1cscs.coefficient({C: 2})
cB = theta1cscs.coefficient({C: 1})
cC = theta1cscs.coefficient({C: 0})
Delta = cB*cB-4*cA*cC
Delta1 = 2*(t_v**3*t_w + 2*t_v**2*t_w**2 - t_v**2*t_w + t_v*t_w - 2*t_v - t_w)
Delta1carre = Delta1 * Delta1
Delta2 = Delta // Delta1carre
print("check: ", Delta - Delta2 * Delta1carre)
check_zero(Delta - Delta2 * Delta1carre)
print("Delta2 = delta20 + delta22 c^2 + delta24 c^4 with")
delta20 = Delta2.coefficient({c: 0}).factor()
delta22 = Delta2.coefficient({c: 2}).factor()
delta24 = Delta2.coefficient({c: 4}).factor()
print("delta20 = ")
print(delta20)
print("delta22 = ")
print(delta22)
print("delta24 = ")
print(delta24)
print("check: ", Delta2 - delta20 - c**2*delta22 - c**4*delta24)
check_zero(Delta2 - delta20 - c**2*delta22 - c**4*delta24)
print("We have Delta2(c^2) >= min( delta20, delta20 + delta22 + delta24, delta20 + delta22/2 )")

print("delta20 = ", delta20, ">= 0")
dummy2 = delta20 + delta22 + delta24
print("delta20 + delta22 + delta24 = ", dummy2.factor(), ">= 0")
dummy3 = 2*delta20 + delta22
print("2*delta20 + delta22 = ")
print(dummy3.factor())
print(" = (1-t_w)*(1-t_v**2)*(1-t_v*t_w)**2 * Lambda(t_v,t_w) with ")
Lambda = dummy3 // ((1-t_w)*(1-t_v**2)*(1-t_v*t_w)**2)
print(" Lambda = t_v^2*(", Lambda.coefficient({t_v: 2}).factor(), \
      "+ t_v*", Lambda.coefficient({t_v: 1}).factor(), \
      "+ (", Lambda.coefficient({t_v: 0}).factor(), ")")

print("We prove that the polynomial Lambda is non negative "
      "under the assumptions -1<t_v<1 and -1<t_w<t_w^*(t_v). "
      "We compute the discriminant of Lambda")
DeltaLambda = Lambda.coefficient({t_v: 1})*Lambda.coefficient({t_v: 1}) - 4*Lambda.coefficient({t_v: 0})*Lambda.coefficient({t_v: 2})
print("Delta_Lambda = ", DeltaLambda.factor())

print("If 0 <= t_w < 1, Delta_Lambda <= 0 implies Lambda >= 0 as ")
print("Lambda(t_v=0) = ", Lambda.subs({t_v: 0}), " > 0.")

print("If -1 < t_w < 0, Delta_Lambda > 0 and Lambda has two distinct real roots.")
print("Lambda(t_v=-1) = ", Lambda.subs({t_v: -1}).factor(), " > 0")
print("Lambda(t_v=+1) = ", Lambda.subs({t_v: 1}).factor(), " < 0")

print("We write t_v^square(t_w) for the zero of Lambda that lies in (-1,1) whenever -1 < t_w < 0. "
      "Thus, it suffices to show that "
      "Lambda(t_v^*(t_w),t_w) >= 0.")
Upsilon = -t_v**2*t_w**2 + 4*t_v**2*t_w + 6*t_v*t_w**2 + t_v**2 - t_w**2 - 6*t_v - 4*t_w + 1
a_x = Upsilon.coefficient({t_v: 2})
b_x = Upsilon.coefficient({t_v: 1})
c_x = Upsilon.coefficient({t_v: 0})
delta_x = b_x*b_x - 4*a_x*c_x
tvstar = 2*c_x / (sa.sqrt(delta_x) - b_x)
print("t_v^*(t_w) = ", tvstar)
racine = sa.var("__rr__")  # symbolic variable for sqrt(8+8*t_w^4)
print("Define __rr__ the symbolic square root of 8(1+t_w^4)")
print("tvstar = (1-4*t_w-t_w^2) / (3*(1-t_w^2) + racine)")
numstar = 1-4*t_w-t_w**2
denomstar = 3*(1-t_w**2) + racine
print("Compute the sign of Lambda(tvstar(tw), tw)")
print("Lambda(tvstar(tw), tw) * (3*(1-t_w**2) + racine)^2 = ")
Lstar = Lambda.coefficient({t_v: 2}) * numstar**2 + Lambda.coefficient({t_v: 1}) * numstar * denomstar + Lambda.coefficient({t_v: 0}) * denomstar**2
Lstar = Lstar.expand()
print(Lstar)
Lstar2 = Lstar.subs(racine**2 == 8+8*t_w**4)
print("check: ", (Lstar-Lstar2).factor())
print(Lstar2.factor())

terme1 = 3*t_w**4  + 4*t_w**3  - 2*t_w**2  - 4*t_w + 3
terme2 = (1 - t_w - t_w**2)
print("3t_w^4+4t_w^3-2t_w^2-4t_w+3 = 1+t_^4+2(1-t_w-t_w^2)^2 > 0")
print("check: ", terme1 - (1+t_w**4 + 2*(1-t_w-t_w**2)**2))
check_zero(terme1 - (1+t_w**4 + 2*(1-t_w-t_w**2)**2))
terme1carremoinsterme2carre = terme1*terme1 - 8*(1+t_w**4)*terme2*terme2
print("(3t_w^4+4t_w^3-2t_w^2-4t_w+3)^2-8(1+t_w^4)(1-t_v-t_w^2)^2 = ")
print(terme1carremoinsterme2carre.factor(), " >= 0.")

end_of_verification()