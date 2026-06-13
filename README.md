# Computer-assisted verification for A new class of Euler explosions

This repository accompanies the paper *A new class of Euler explosions*
([arXiv](https://www.arxiv.com)) and contains the code used to reproduce the
computer-assisted parts of the proof.

## Requirements

The repository contains two independent verification components.

The rational-cover verification must be run in a Python environment with
SymPy available. In particular, the following import must work:

    import sympy

The interval-bound verification must be run in a SageMath-enabled Python
environment. In particular, the following import must work:

    from sage.all import *

## Rational cover

All files for this part are located in the `rational_cover/` directory.

### `verify_bernstein_q_rational_cover.py`

Exact rational-cover verification for `lemma:Bernstein:q`.

  * Starts from the phase-plane definitions of `P_v`, `P_q`, and the quadratic
    barrier.
  * Derives the separately affine forms in
    `scaled_c = (1+3*alpha)c_r`, `scaled_v = (1+3*alpha)V_2`, and `rr < 0`.
  * Checks the small-alpha case.
  * Checks the rational cover of `[10^(-2),1+sqrt(2)]` using exact rational
    arithmetic.
  * Uses SymPy only to take exact floors and ceilings of the algebraic endpoint
    values defining the rational boxes.

Running this script reproduces the rational-cover verification:

    cd rational_cover
    python verify_bernstein_q_rational_cover.py

Optional flags:

  * `--box-digits N` changes the number of decimal digits used in the rational
    endpoint boxes.
  * `--bernstein-degree N` uses one common Bernstein degree for interval
    polynomials.
  * `--tex-tables tables.tex` creates or overwrites that file with the
    generated TeX tables. If a path includes a directory, the directory must
    already exist.
  * `--verbose` prints the normal-form details and every interval in the cover.

## Computer-assisted bounds

All files for this part are located in the `computer_assisted_bounds/`
directory.

### `script.py`

SageMath interval verifier for the bounds displayed in
`Table~\ref{table:bb}`.

The main entry point is `script.py`. Run it from a SageMath-enabled Python
environment. The following import must work:

    from sage.all import *

The script verifies the six two-sided interval bounds for
`\tilde{\BB}_0,\ldots,\tilde{\BB}_5` reported in
`Table~\ref{table:bb}`. The verification is performed with Sage's
`RealBallField` and `RealIntervalField` arithmetic.

The variable used internally is

    b = sqrt(alpha),     b in [0, 14/9].

For each of the six rows, the verifier:

  * builds the exact radical expression from the formulas in `lemmas.py`;
  * reads the claimed lower and upper table bounds from
    `supplementary_data/bounds.txt`;
  * subdivides the interval using the cuts configured in `parameters.py`;
  * proves that every resulting subinterval lies strictly inside the claimed
    lower and upper bounds;
  * reports failure if any subinterval cannot be certified before the
    configured depth or box limit.

The printed decimal traces are only diagnostics. The mathematical check is the
interval comparison made in `lemmas.py`.

After activating any Python environment in which `from sage.all import *`
works, run:

    cd computer_assisted_bounds
    python script.py

A successful run ends with:

    The bounds reported in Table~\ref{table:bb} have been verified.

Useful flags:

  * `--verbose` or `-v` prints row-level progress and the bounds read from
    `bounds.txt`.
  * `--print-subintervals` prints every accepted subinterval and its local
    margins.
  * `-vv` is equivalent to verbose mode with the subinterval trace enabled.

Files in this component:

  * `script.py` is the command-line driver. It runs the table-bound
    verification and the consistency check for the bound labels.
  * `lemmas.py` is the interval verifier. It defines the expression DSL, the
    explicit formulas for `b5b0` through `b5b5`, the adaptive subdivision
    routine, and the lemma-level checks.
  * `parameters.py` contains the global numerical settings: 200-bit ball
    precision, the interval `[0, 14/9]`, initial cuts, maximum subdivision
    depth, box limit, and printing controls.
  * `supplementary_data/bounds.txt` contains the rounded lower and upper
    bounds appearing in `Table~\ref{table:bb}`. Each non-comment line has the
    form `label|value|side|description`.

If a table bound changes, update `supplementary_data/bounds.txt` and rerun
`script.py`. The labels and their `lower`/`upper` side tags are checked by
`substituting_estimates()`, so a missing label or a swapped side is reported as
a verification failure.
