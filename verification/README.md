# Symbolic verification of stability results

This directory contains SageMath scripts implementing symbolic calculations associated with the stability analysis of a D1Q3 lattice Boltzmann scheme for the advection equation.

The scripts aim to make the algebraic calculations underlying the mathematical results reproducible. They use polynomial manipulations, Cohn–Schur transformations, and calculations involving the amplification matrix.

## Requirements

The scripts require SageMath. They are intended to be executed using SageMath's Python interpreter.

Run the commands below from the root directory of the repository.

## Scripts

### `utils.py`

This module contains shared utilities used by the verification scripts, including:

* construction and manipulation of polynomial expressions;
* Cohn–Schur transformations;
* polynomial differentiation and substitution;
* construction of the amplification matrix;
* common functions for reporting verification results.

### `expressions_of_the_polynomials.py`

This script checks explicit polynomial expressions arising from the algebraic transformations used in the stability analysis. It compares expressions obtained symbolically with their expected forms.

### `verification_section_5_2.py`

This script implements symbolic calculations associated with Section 5.2 of the paper, including Lemmas 20–23. The calculations use polynomial transformations and algebraic identities to study the stability conditions.

### `verification_prop27.py`

This script implements symbolic calculations associated with Proposition 27. It studies stability conditions in the interior of the parameter domain using polynomial transformations and algebraic sign conditions.

### `verification_section_5_4.py`

This script implements calculations associated with the edge-stability results in Section 5.4, including Propositions 28–31. It considers cases in which one of the parameters lies on the boundary of its admissible interval.

### `verification_section_5_5.py`

This script implements calculations associated with the vertex-stability results in Section 5.5, including Propositions 32–34. It studies cases in which both parameters lie on the boundary of their admissible intervals.

### `verification_prop34.py`

This script focuses on an algebraic verification associated with Proposition 34. In particular, it examines the amplification matrix, candidate eigenvectors, and the corresponding algebraic relations in a vertex-stability case.

## Running the scripts

Run an individual script from the root of the repository. For example:

```bash
sage -python verification/expressions_of_the_polynomials.py
```

To run all the verification scripts in sequence:

```bash
for f in \
    verification/expressions_of_the_polynomials.py \
    verification/verification_prop27.py \
    verification/verification_prop34.py \
    verification/verification_section_5_2.py \
    verification/verification_section_5_4.py \
    verification/verification_section_5_5.py
do
    echo "===== $f ====="
    sage -python "$f" || exit 1
done
```

The loop stops if a script exits with a nonzero status. This is useful for detecting execution errors, provided that failed symbolic checks also cause the corresponding script to exit with a nonzero status.

## Interpreting the results

The scripts perform symbolic calculations and check algebraic identities associated with the stability analysis. Depending on the calculation, these checks may involve polynomial coefficients, factorizations, reductions modulo polynomial ideals, or identities involving the amplification matrix.

A successful run indicates that the checks implemented in the script have passed. It should not be interpreted as an independent proof of every mathematical result in the paper: the validity of the conclusions also depends on the mathematical assumptions and on the correctness of the implementation.

For identities reduced modulo a polynomial ideal, the result establishes the identity in the corresponding quotient ring, subject to the relations defining that ideal.

## Relationship with the paper

The scripts are organized around the polynomial calculations and the stability cases studied in the paper. The names of the scripts indicate the relevant sections or propositions.

For the full scientific context, consult the manuscript in [`../paper/manuscrit.pdf`](../paper/manuscrit.pdf) and the main repository [`../README.md`](../README.md).
