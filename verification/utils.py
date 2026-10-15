"""
Utility functions for manipulating 
polynomial coefficients and amplification matrices.

This module gathers a few helper routines 
used in the SageMath computations related to
complex polynomial coefficients, symbolic simplification, 
and the characteristic
polynomial of a 3x3 amplification matrix.
"""

import sage.all as sa  # type: ignore[import-untyped]


def banner(message: str) -> None:
    """Display a message in a framed banner."""
    lines = message.splitlines()
    width = 75  #max(len(line) for line in lines)

    border = "=" * (width + 4)

    print()
    print(border)
    for line in lines:
        print(f"  {line.ljust(width)}  ")
    print(border)
    print()


def check_zero(*args) -> None:
    """Check that symbolic expressions are identically zero."""
    for expr in args:
        if expr != 0:
            print("ERROR: non zero residual")
            print(f"Residual: {expr}")
            raise AssertionError("The expression is not zero.")


def end_of_verification() -> None:
    """Display a message indicating that the verification completed successfully."""
    message = "Verification completed successfully."
    banner(message)


def CohnSchur(Pre, Pim):
    """Compute the Cohn-Schur relations for a pair of coefficient arrays.

    Args:
        Pre: Real parts of the polynomial coefficients, ordered from degree 0 to d.
        Pim: Imaginary parts of the polynomial coefficients, ordered from degree 0 to d.

    Returns:
        A pair (Qre, Qim) containing the real and imaginary parts of the derived
        polynomial coefficients.
    """
    Qre, Qim = [], []
    d = len(Pre) - 1
    for k in range(1, d+1):
        dummy = Pre[d]*Pre[k] + Pim[d]*Pim[k] - Pre[0]*Pre[d-k] - Pim[0]*Pim[d-k]
        Qre.append(dummy)
        dummy = Pre[d]*Pim[k] - Pim[d]*Pre[k] + Pre[0]*Pim[d-k] - Pim[0]*Pre[d-k]
        Qim.append(dummy)
    return Qre, Qim


def Derive(Pre, Pim):
    """Differentiate a polynomial represented by its real and imaginary coefficients.

    The coefficient lists are assumed to be ordered from degree 0 to the maximum
    degree. The derivative is computed coefficient-wise as k * a_k for each degree k.

    Args:
        Pre: Real parts of the polynomial coefficients.
        Pim: Imaginary parts of the polynomial coefficients.

    Returns:
        The derivative coefficients as a pair (Qre, Qim).
    """
    Qre, Qim = [], []
    d = len(Pre) - 1
    for k in range(1, d+1):
        Qre.append(k*Pre[k])
        Qim.append(k*Pim[k])
    return Qre, Qim


def simp_substitute(expr, exprin, exprout):
    """Simplify an expression by repeatedly removing occurrences of ``exprin``.

    This routine performs a basic substitution-style simplification using integer
    division in a loop to eliminate factors of ``exprin`` and replace them by
    ``exprout``. It is mainly used for symbolic cleanup in polynomial identities.

    Args:
        expr: The expression to simplify.
        exprin: The sub-expression to remove.
        exprout: The replacement expression.

    Returns:
        The simplified expression.
    """
    # return expr.substitute(
    #     {
    #         exprin: exprout
    #     }
    # )
    coeff = expr // exprin
    while coeff != 0:
        expr = expr - (exprin-exprout) * coeff
        coeff = expr // exprin
    return expr


def extract_re_im(expr, reii, imii, degre):
    """Separate a symbolic expression into its real and imaginary parts.

    Args:
        expr: The expression to decompose.
        reii: Symbol representing the real part parameter.
        imii: Symbol representing the imaginary part parameter.
        degre: Maximum polynomial degree used in the decomposition.

    Returns:
        A tuple (expr_re, expr_im) with the real and imaginary contributions.
    """
    expr_re = expr.subs(reii=1, imii=0)
    expr_im = 0
    for n in range(1, degre+1):
        p, r = n // 2, n % 2
        coeff = expr.coefficient({imii: n}).subs(reii=1)
        if r == 0:
            expr_re += (-1)**p * coeff
        else:
            expr_im += (-1)**p * coeff
    return expr_re, expr_im


def amplification_matrix(c, tv, tw, C, S, reii, imii, print_char=False):
    """Build the 3x3 amplification matrix associated with the model parameters.

    Args:
        c: Courant-like parameter.
        tv: First time-variation parameter.
        tw: Second time-variation parameter.
        C: Cosine-related quantity used in the symbolic model.
        S: Sine-related quantity used in the symbolic model.
        reii: Symbol for the real part.
        imii: Symbol for the imaginary part.
        print_char: If True, print the characteristic polynomial and its real/imaginary
            coefficients for inspection.

    Returns:
        The SageMath matrix describing the amplification process.
    """
    sv, sw = 1+tv, 1+tw
    d = c*c
    G = sa.matrix(
        [[
           reii*(1-sw*d),  reii*sw*(1-d), reii*sw*(1-d)
        ], [
            (C*reii-S*imii)/2*(sw*d+sv*c),
            (C*reii-S*imii)/2*(2-sw*(1-d)-sv*(1-c)),
            (C*reii-S*imii)/2*(-sw*(1-d)+sv*(1+c))
        ], [
            (C*reii+S*imii)/2*(sw*d-sv*c),
            (C*reii+S*imii)/2*(-sw*(1-d)+sv*(1-c)),
            (C*reii+S*imii)/2*(2-sw*(1-d)-sv*(1+c))
        ]]
    )

    if print_char:
        XX = sa.var('X', domain=sa.RR)
        P = G.charpoly(XX)
        liste_coeff = P.coefficients(sparse=False)
        liste_coeff_re, liste_coeff_im = [], []
        for k in range(len(liste_coeff)):
            dummy = liste_coeff[k]
            cre, cim = extract_re_im(dummy, reii, imii, len(liste_coeff)-1)
            cre = simp_substitute(cre, 1-S**2, C**2)
            cim = simp_substitute(cim, 1-S**2, C**2)
            liste_coeff_re.append(cre)
            liste_coeff_im.append(cim)
        print("*"*80)
        print("Characteristic polynomial")
        for k in range(len(liste_coeff)):
            print("-"*80)
            print(f"Degree {k}")
            print("Real part")
            print(liste_coeff_re[k])
            print("Imaginary part")
            print(liste_coeff_im[k])
        print("*"*80)

    return G
