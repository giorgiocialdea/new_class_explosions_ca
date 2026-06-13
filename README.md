# Computer-assisted verification for A new class of Euler explosions

This repository accompanies the paper [A new class of Euler explosions](https://www.arxiv.com) and contains the code and supplementary data used to reproduce the computer-assisted parts of the proof.

## Requirements

The rational-cover verification must be run in a Python environment with SymPy. In particular, the following import must work:

    import sympy

The interval-arithmetic verification must be run in a SageMath-enabled Python environment. In particular, the following import must work:

    from sage.all import *

## Python source files

### `rational_cover/verify_bernstein_q_rational_cover.py`

Script verifying the rational-cover part of Lemma~\ref{lemma:Bernstein:q}.

  * Uses exact rational arithmetic.
  * Checks the small-alpha case and the rational cover of `[10^(-2),1+sqrt(2)]`.
  * Can optionally print the generated TeX tables.

Running this script reproduces the rational-cover verification:

    cd rational_cover
    python verify_bernstein_q_rational_cover.py

Optional flags:

  * `--box-digits N` changes the number of decimal digits used in the rational endpoint boxes.
  * `--bernstein-degree N` uses one common Bernstein degree for the interval polynomials.
  * `--tex-tables tables.tex` creates or overwrites the named file with the generated TeX tables. If the path includes a directory, that directory must already exist.
  * `--verbose` prints the normal-form details and every interval in the cover.

### `computer_assisted_bounds/script.py`

Main script reproducing the interval-arithmetic verification of the bounds reported in Table~\ref{table:bb}.

  * Calls the verification routines defined in `lemmas.py`.
  * Outputs timing information and verification status.

After activating any Python environment in which `from sage.all import *` works, run:

    cd computer_assisted_bounds
    python script.py

A successful run ends with:

    The bounds reported in Table~\ref{table:bb} have been verified.

Optional flags:

  * `--verbose` or `-v` prints row-level progress and the bounds read from `bounds.txt`.
  * `--print-subintervals` prints every accepted subinterval and its local margins.
  * `-vv` is equivalent to verbose mode with the subinterval trace enabled.

### `computer_assisted_bounds/lemmas.py`

Verification routines for the interval-arithmetic bounds in Table~\ref{table:bb}.

### `computer_assisted_bounds/parameters.py`

Defines global numerical parameters used in the interval-arithmetic computations.

  * `VERBOSE` controls the amount of printed information.

## Data files

### `computer_assisted_bounds/supplementary_data/bounds.txt`

  * Read with `load_bounds()`.
  * Contains the bounds reported in Table~\ref{table:bb}.
