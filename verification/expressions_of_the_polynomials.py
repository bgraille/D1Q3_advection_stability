r"""
Expression of the polynomials

Parameters
----------
- c   Courant-like parameter
- t_v first rescaled relaxation parameter $t_v=s_v-1$;
- t_w second rescaled relaxation parameter $t_w=s_w-1$;
- C   Cosine-related quantity used in the symbolic model $C=\cos(\xi)$; 
- S   Sine-related quantity used in the symbolic model $S=\sin(\xi)$.
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


def verification(Pre, Pim, Pre_rw, Pim_rw):
    """check if the two polynomials are identical"""
    for k in range(len(Pre)):
        print(f"Degree {k}")
        print("-"*8)
        print("Real part")
        try:
            Pre[k] = simp_substitute(Pre[k], S**2, 1-C**2)
            print(Pre[k].factor())
        except (ArithmeticError, AttributeError, TypeError):
            print(Pre[k])
        print("Verification: ", Pre[k]-Pre_rw[k])
        check_zero(Pre[k]-Pre_rw[k])
        print("Imaginary part")
        try:
            print(Pim[k].factor())
        except (ArithmeticError, AttributeError, TypeError):
            print(Pim[k])
        print("Verification: ", Pim[k]-Pim_rw[k])
        check_zero(Pim[k]-Pim_rw[k])
        print()


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
# Compute the polynomials
##############################################################################
P2cs_re, P2cs_im = CohnSchur(P3_re, P3_im)
P2d_re, P2d_im = Derive(P3_re, P3_im)
P1cscs_re, P1cscs_im = CohnSchur(P2cs_re, P2cs_im)
P1csd_re, P1csd_im = Derive(P2cs_re, P2cs_im)
P1dcs_re, P1dcs_im = CohnSchur(P2d_re, P2d_im)
P1dd_re, P1dd_im = Derive(P2d_re, P2d_im)


##############################################################################
# P2cs
##############################################################################
banner("Polynomial P2cs after one step of Cohn-Schur algorithm")
P2cs_re_rw = [
    -c**2*t_v*(1-t_w**2)*(1-C)-C*(t_v+t_w)*(1-t_v*t_w),
    c**2*(1+t_w)*(1-C)*(1-t_v**2*t_w) - (1-t_v**2*t_w**2) + C*t_v*(1-t_w**2) + C*t_w*(1-t_v**2)
,
    1-t_v**2*t_w**2
]
P2cs_im_rw = [
    (1-t_v**2)*c*S*t_w,
    c*(1-t_v*t_w**2)*(1+t_v)*S,
    0
]
verification(
    P2cs_re, P2cs_im,
    P2cs_re_rw, P2cs_im_rw
)


##############################################################################
# P2d
##############################################################################
banner("Polynomial P2d after one step of derivation")
P2d_re_rw = [
    t_v*t_w*(1-c**2*(1-C)) - t_v*((1-C)*c**2+C) - t_w*C,
    -2*(1-c**2*(1-C)) + 2*t_w*((1-C)*c**2+C) + 2*t_v*C,
    3
]
P2d_im_rw = [
    c*t_w*(1+t_v)*S,
    2*c*(1+t_v)*S,
    0
]
verification(
    P2d_re, P2d_im,
    P2d_re_rw, P2d_im_rw
)


##############################################################################
# P1cscs
##############################################################################
banner("Polynomial P1cscs after two steps of Cohn-Schur algorithm")
P1cscs_re_rw = [
    - (1-t_v*t_w)**2 * ((1+t_v*t_w)**2-C**2*(t_v+t_w)**2)
    + c**2*(1-C)*(
        (1-t_v)*(
            (1 - t_v*t_w)*(1+t_v**2*t_w**3) + t_v*t_w*(t_w-t_v)*(1+t_w)
        ) + C*(1+t_v) * (
            (t_v+t_w**2)*(1-t_v*t_w)**2 + t_v*t_w*(1-t_w)**2*(1+t_v)
        )
    ) + c**4*(1-t_w)*t_v*(1-C)**2*(1+t_w)**2*(1-t_v**2*t_w),
    (1-t_v*t_w)**2*((1+t_v*t_w)**2-C**2*(t_v+t_w)**2)
    + (1-C)*c**2*(
        -(1+C)*t_w**2*(1-t_v**2)**2
        - 2*C*(1-t_w**2)*t_v*(t_v+t_w)*(1-t_v*t_w)
    )
    -c**4 * (C - 1)**2 * (t_w - 1)**2 * (t_w + 1)**2 * t_v**2
]
P1cscs_im_rw = [
    S * (t_w + 1) * (t_v + 1) * c * (
        (1-t_v*t_w)**2*(1+t_v*t_w-C*t_v-C*t_w)
        - c**2*(1-C)*(
            t_v*(1-t_w)*(1-t_v*t_w**2)
            + t_w*(1-t_v)*(1-t_w*t_v**2)
        )
    ),
    0
]
verification(
    P1cscs_re, P1cscs_im,
    P1cscs_re_rw, P1cscs_im_rw
)


##############################################################################
# P1csd
##############################################################################
banner(
    "Polynomial P1csd after\n"
    "- one step of Cohn-Schur algorithm \n"
    "- one step of derivation"
)
P1csd_re_rw = [
    c**2*(1+t_w)*(1-C)*(1-t_v**2*t_w) - (1-t_v**2*t_w**2) + C*t_v*(1-t_w**2) + C*t_w*(1-t_v**2)
,
    2*(1-t_v**2*t_w**2)
]
P1csd_im_rw = [
    c*(1-t_v*t_w**2)*(1+t_v)*S,
    0
]
verification(
    P1csd_re, P1csd_im,
    P1csd_re_rw, P1csd_im_rw
)


##############################################################################
# P1dcs
##############################################################################
banner(
    "Polynomial P1dcs after\n"
    "- one step of derivation\n"
    "- one step of Cohn-Schur algorithm"
)
P1dcs_re_rw = [
    -2*(1-C*(t_v+t_w))*(3-t_v*t_w + C*(t_v+t_w)) + 2*(1-C)*c**2*(
        3*(1+t_w) + C*(t_v+t_w**2)*(1+t_v) -t_v*(1+t_w)**2 - t_w*(1+t_v)**2
    ) + 2*c**4*t_v*(1-C)**2*(1+t_w)**2,
    9-(C*t_v+C*t_w-t_v*t_w)**2 - c**4*(1-C)**2*(1+t_w)**2*t_v**2 + (1-C)*c**2*(
        2*t_v*(1+t_w)*(t_w*t_v-C*t_v-C*t_w) - (1+C)*t_w**2*(1+t_v)**2
    )
]
P1dcs_im_rw = [
    2*c*S*(1+t_v)*(
        3+t_w*(1+t_v)-(C+c**2*(1-C))*(t_v+t_w)*(1+t_w)
    ),
    0
]
verification(
    P1dcs_re, P1dcs_im,
    P1dcs_re_rw, P1dcs_im_rw
)


##############################################################################
# P1dd
##############################################################################
banner("Polynomial P1dd after two steps of derivation")
P1dd_re_rw = [
    -2*(1-c**2*(1-C)) + 2*t_w*((1-C)*c**2+C) + 2*t_v*C,
    6
]
P1dd_im_rw = [
    2*c*(1+t_v)*S,
    0
]
verification(
    P1dd_re, P1dd_im,
    P1dd_re_rw, P1dd_im_rw
)


end_of_verification()
