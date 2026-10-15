# Computer-assisted verification of stability results for a D1Q3 Lattice Boltzmann scheme

This repository contains the computer-assisted calculations associated with the paper *Numerical analysis of a Lax-Wendroff type lattice Boltzmann scheme for the advection equation*.

The project provides SageMath scripts for the symbolic verification of mathematical results and a Python animation for the interactive exploration of the eigenvalues of the amplification matrix.

## Repository contents

```text
.
├── CITATION.cff
├── LICENSE
├── README.md
├── pyproject.toml
├── paper/
│   └── manuscrit.pdf
├── scripts/
│   └── check_environment.py
├── verification/
│   ├── README.md
│   ├── expressions_of_the_polynomials.py
│   ├── utils.py
│   ├── verification_prop27.py
│   ├── verification_prop34.py
│   ├── verification_section_5_2.py
│   ├── verification_section_5_4.py
│   └── verification_section_5_5.py
└── visualization/
    └── eigenvalues_visualization.py
```

* `paper/`: the manuscript associated with the verification scripts.
* `verification/`: SageMath scripts implementing symbolic calculations related to the results in the paper.
* `visualization/`: an interactive animation for exploring the eigenvalues of the amplification matrix.
* `scripts/`: utilities for checking the computational environment.
* `CITATION.cff`: citation metadata.
* `LICENSE`: the project license.

## Requirements

The symbolic verification scripts require SageMath.

The interactive visualization requires Python, NumPy, and Matplotlib. These packages must be available in the Python environment used to run the animation.

## Checking the computational environment

Run the following command from the root of the repository:

```bash
sage -python scripts/check_environment.py
```

The script reports the versions of Python, SageMath, NumPy, and Matplotlib available in the environment and checks that the required modules can be imported.

## Symbolic verification

The scripts in `verification/` perform symbolic calculations associated with the stability analysis of the D1Q3 lattice Boltzmann scheme.

The main scripts are:

| Script                              | Purpose                                                                                                      |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `expressions_of_the_polynomials.py` | Checks the expressions of polynomials obtained through the algebraic transformations used in the analysis.   |
| `verification_prop27.py`            | Symbolic calculations associated with Proposition 27.                                                        |
| `verification_prop34.py`            | Algebraic verification involving the amplification matrix and its eigenstructure in a vertex-stability case. |
| `verification_section_5_2.py`       | Calculations associated with the results in Section 5.2, including Lemmas 20–23.                             |
| `verification_section_5_4.py`       | Calculations associated with the edge-stability results in Section 5.4, including Propositions 28–31.        |
| `verification_section_5_5.py`       | Calculations associated with the vertex-stability results in Section 5.5, including Propositions 32–34.      |

Run an individual script from the repository root using SageMath. For example:

```bash
sage -python verification/verification_prop27.py
```

To execute all the verification scripts in sequence:

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

The scripts use the shared utilities defined in `verification/utils.py`.

The implemented checks test symbolic identities and algebraic relations arising in the stability analysis. A successful execution means that the implemented checks have passed; it does not, by itself, establish the correctness of every mathematical statement in the paper independently of the assumptions and implementation.

Further details on the verification scripts are provided in [`verification/README.md`](verification/README.md).

## Interactive visualization

The script `visualization/eigenvalues_visualization.py` provides an interactive animation for exploring the eigenvalues of the amplification matrix as parameters vary.

The user can manipulate the slider to investigate parameter values and observe the relationship between the eigenvalues and the stability condition.

Run the animation from the repository root:

```bash
python visualization/eigenvalues_visualization.py
```

Use a Python environment in which the dependencies required by the script are installed.

The animation is intended for numerical exploration and illustration. It does not replace the symbolic verifications or constitute a proof of stability.

## Scientific reference

The scripts accompany the paper:

*Numerical analysis of a Lax-Wendroff type lattice Boltzmann scheme for the advection equation.*

See [`paper/manuscrit.pdf`](paper/manuscrit.pdf) for the manuscript included in this repository.

## Authors

* Thomas Bellotti
* Filipa Caetano
* François Dubois
* Benjamin Graille

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Please use the bibliographic information provided there when citing this software.

## License

This project is distributed under the MIT License. See [`LICENSE`](LICENSE) for details.
