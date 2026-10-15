"""Check the Python environment used by the project.

This script reports the versions of Python, SageMath, NumPy, and Matplotlib
available in the current environment.
"""

import sys

import matplotlib
import numpy
import sage.all as sage  # type: ignore[import-untyped]


def main() -> None:
    """Display the versions of the required software."""
    print("-" * 80)
    print("Computational environment")
    print("-------------------------")
    print(f"Python:     {sys.version.split()[0]}")
    print(f"SageMath:   {sage.version()}")
    print(f"NumPy:      {numpy.__version__}")
    print(f"Matplotlib: {matplotlib.__version__}")
    print()
    print("Environment check completed successfully.")
    print("-" * 80)


if __name__ == "__main__":
    main()
