r"""
Verification script for the section 5.2

Parameters
----------
- c   Courant-like parameter
- t_v first rescaled relaxation parameter t_v = s_v - 1;
- t_w second rescaled relaxation parameter t_w = s_w - 1;
- C   Cosine-related quantity used in the symbolic model C = cos(xi); 
- S   Sine-related quantity used in the symbolic model S = sin(xi).

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
# __Proof of Lemma 20__
##############################################################################
banner("Proof of Lemma 20")
dummy = (P3_re[3]**2+P3_im[3]**2) - (P3_re[0]**2+P3_im[0]**2)
print(r"|\alpha_{3,3}|^2 - |\alpha_{3,0}|^2")
print(" = ")
print(dummy)
print(" = ")
print(dummy.factor(), " > 0")

##############################################################################
# __Proof of Lemma 21__
##############################################################################
banner("Proof of Lemma 21")
print("Compute P2cs - polynomial after one step of Cohn-Schur algorithm")
P2cs_re, P2cs_im = CohnSchur(P3_re, P3_im)
print("P2cs = 0 <=> ")
for eq in P2cs_re + P2cs_im:
    if eq != 0:
        print("   ", eq.factor(), " = 0")
print(r"The condition \alpha_{2,2}^\circ = 0 reads")
print(P2cs_re[2], " = 0 <=> |t_v t_w| = 1 <=> |t_v| = |t_w| = 1")
print(r"check all combination t_v, t_w \in {-1,1}")
for valtv in [-1, 1]:
    for valtw in [-1, 1]:
        dummy = [ck.subs({t_v: valtv, t_w: valtw}) for ck in P2cs_re]
        print("check the real part ", *dummy)
        check_zero(*dummy)
        dummy = [ck.subs({t_v: valtv, t_w: valtw}) for ck in P2cs_im]
        print("check the imaginary part ", *dummy)
        check_zero(*dummy)
print("Then |t_v| = |t_w| = 1 <=> P2cs = 0")

##############################################################################
# __Proof of Lemma 22__
##############################################################################
banner("Proof of Lemma 22")
print("Compute P2cs - polynomial after one step of Cohn-Schur algorithm")
P2cs_re, P2cs_im = CohnSchur(P3_re, P3_im)
print("P2cs(0) = (", P2cs_re[0], ", ", P2cs_im[0], ")")
print("P2cs*(0) = (", P2cs_re[2], ", ", -P2cs_im[2], ")")
print("\n|P2cs(0)| < |P2cs*(0)| <=> cond > 0 \nwhere cond = ")
cond = P2cs_re[2]**2 - P2cs_re[0]**2-P2cs_im[0]**2
cond = simp_substitute(cond, S**2, 1-C**2)
print(cond)
print("cond is a quadratic polynomial in c^2 that reads ")
alpha1 = []
for k in range(cond.degree(c)+1):
    dummy = cond.coefficient({c: k})
    alpha1.append(dummy)
for k, alpha1k in enumerate(alpha1):
    if alpha1k != 0:
        print(f"Terme of degree {k}")
        print(alpha1k.factor())
print("""
The coefficient of c^4 is not positive
cond is then a concave polynomial in c^2.
""")

print("First, we consider |t_v|<1 and |t_w|<1")
print("Evaluation in c^2=0")
print(alpha1[0].factor())
print(" = (1-t_v*t_w)**2*((1-t_v**2)*(1-t_w**2)+(1-C**2)*(t_v+t_w)**2)")
alpha10 = (1-t_v*t_w)**2*((1-t_v**2)*(1-t_w**2)+(1-C**2)*(t_v+t_w)**2)
print("Verification: ", alpha1[0] - alpha10)
print("As |t_v|<1 and |t_w|<1, this term is positive!")
print("Evaluation in c^2=1")
dummy = alpha1[0] + alpha1[2] + alpha1[4]
print(dummy.factor())
print(" = (1-t_v**2)*(1-t_w**2) * ( (1-C*t_v*t_w)**2 + (1-C**2)*t_v**2*t_w**2 )")
dummy2 = (1-C*t_v*t_w)**2 + (1-C**2)*t_v**2*t_w**2
dummy_new = dummy2 * (1-t_v**2)*(1-t_w**2)
print("Verification: ", dummy - dummy_new)
print("As |t_v|<1 and |t_w|<1, this term is positive!")

print("Second, we consider |t_v|=1 and |t_w|<1")
# tv = -1
cond_edgestt_l = simp_substitute(cond, t_v, -1)
print(cond_edgestt_l.factor())
# tv = 1
cond_edgestt_r = simp_substitute(cond, t_v, 1)
print(cond_edgestt_r.factor())
print("> 0 <=> C<1 and c<1 and (c>0 or C>-1)")

print("Third, we consider |t_v|<1 and |t_w|=1")
# tw = -1
cond_edgestt_b = simp_substitute(cond, t_w, -1)
print(cond_edgestt_b.factor())
# tw = 1
cond_edgestt_t = simp_substitute(cond, t_w, 1)
print(cond_edgestt_t.factor())
print("> 0 <=> -1<C<1 and c<1")

##############################################################################
# __Proof of Lemma 23__
##############################################################################
banner("Proof of Lemma 23")
print("Compute P1cscs - polynomial after two steps of Cohn-Schur algorithm")
P1cscs_re, P1cscs_im = CohnSchur(P2cs_re, P2cs_im)

print("First, we consider t_v=-1 and |t_w|<1")
val = -1
print(f"t_v = {val}")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val})
    if terme == 0:
        print(terme)
    else:
        print(terme.factor())
print("All terms cancel out iff C=1 or c=1 or (C=-1 and c=0)")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val, C: 1, S: 0})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val, c: 1})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val, C: -1, S: 0, c: 0})
    check_zero(terme)
print("check: 0")

print("Second, we consider t_v=+1 and |t_w|<1")
val = 1
print(f"t_v = {val}")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val})
    if terme == 0:
        print(terme)
    else:
        print(terme.factor())
print("All terms cancel out iff C=1 or c=1 or (C=-1 and c=0)")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val, C: 1, S: 0})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val, c: 1})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_v: val, C: -1, S: 0, c: 0})
    check_zero(terme)
print("check: 0")

print("Third, we consider |t_v|<1 and t_w=-1")
val = -1
print(f"t_w = {val}")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val})
    if terme == 0:
        print(terme)
    else:
        print(terme.factor())
print("All terms cancel out iff C^2=1 or c=1")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val, C: 1, S: 0})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val, C: -1, S: 0})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val, c: 1})
    check_zero(terme)
print("check: 0")

print("Fourth, we consider |t_v|<1 and t_w=+1")
val = 1
print(f"t_w = {val}")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val})
    if terme == 0:
        print(terme)
    else:
        print(terme.factor())
print("All terms cancel out iff C^2=1 or c=1")
for coeff in [P1cscs_re[0], P1cscs_re[1], P1cscs_im[0], P1cscs_im[1]]:
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val, C: 1, S: 0})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val, C: -1, S: 0})
    check_zero(terme)
    terme = simp_substitute(coeff, S**2, 1-C**2).subs({t_w: val, c: 1})
    check_zero(terme)
print("check: 0")

end_of_verification()
