#!/usr/bin/env python3
"""SymPy-style rational-cover verification for Lemma Bernstein:q.

The script verifies two ranges:

1. The small-alpha range 0 < alpha <= 10^(-2), using x=sqrt(alpha).
2. One unified rational cover from 10^(-2) to 2.4143.  Since 2.4143 is slightly
   larger than 1+sqrt(2), this covers the full interval up to 1+sqrt(2).

The scaled variables are

    C = (1+3 alpha)c_r,
    Y = (1+3 alpha)V_2,
    rr is the paper variable.

In the explosion-side convention rr is negative.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# How this file is organized
# ---------------------------------------------------------------------------
# 1. Declare the symbols and the rational cover of [0,1+sqrt(2)].
# 2. Use SymPy to derive the six Bernstein polynomials L_j.
# 3. Store each L_j in a sparse/dense coefficient format that is quick to
#    evaluate at rational box corners.
# 4. Build rational boxes for C,Y,rr on each cover interval.
# 5. Prove positivity by converting the polynomials in alpha to
#    Bernstein form on each subinterval and checking that all such coefficients
#    have positive lower bounds


import argparse
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from pathlib import Path
from time import perf_counter
from typing import Dict, List, Optional, Sequence, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# Symbols and exact phase-plane definitions
# ---------------------------------------------------------------------------

# alpha is the parameter, and t is the barrier parameter.
alpha, t = sp.symbols("alpha t")

# C,Y,rr are the scaled algebraic variables.  V,Q are the phase-plane variables
# used only while building the expressions from vector-fields.
C, Y, rr = sp.symbols("C Y rr")
V, Q = sp.symbols("V Q")

# The rational cover ends at 2.4143, slightly larger than 1 + sp.sqrt(2).
TRUE_ENDPOINT = 1 + sp.sqrt(2)
RATIONAL_ENDPOINT = sp.Rational(24143, 10000)
assert RATIONAL_ENDPOINT > TRUE_ENDPOINT


@dataclass(frozen=True)
class CoverInterval:
    """One closed interval in the rational alpha-cover."""

    left: str
    right: str

    @property
    def a0(self) -> sp.Rational:
        return sp.Rational(self.left)

    @property
    def a1(self) -> sp.Rational:
        return sp.Rational(self.right)

    def tex_interval(self) -> str:
        return f"[{self.left},{self.right}]"


@dataclass(frozen=True)
class CoverResult:
    """Data recorded for one certified interval check."""

    idx: int
    interval: CoverInterval
    minimum: sp.Rational
    coefficient_j: int
    C_bounds: Tuple[sp.Rational, sp.Rational]
    Y_bounds: Tuple[sp.Rational, sp.Rational]
    rr_bounds: Tuple[sp.Rational, sp.Rational]


# The rational cover of [0, 1+sqrt[2]].
COVER_ENDPOINTS: Tuple[str, ...] = (
    "0.0100", "0.0476", "0.0997", "0.1233", "0.1436", "0.1638",
    "0.1858", "0.2116", "0.2435", "0.2765", "0.3047", "0.3294",
    "0.3517", "0.3726", "0.3919", "0.4091", "0.4247", "0.4394",
    "0.4537", "0.4679", "0.4823", "0.4971", "0.5124", "0.5285",
    "0.5455", "0.5546", "0.5637", "0.5832", "0.6042", "0.6269",
    "0.6516", "0.6785", "0.7078", "0.7398", "0.7748", "0.8131",
    "0.8549", "0.9006", "0.9504", "1.0000", "1.0569", "1.1111",
    "1.1630", "1.2129", "1.2610", "1.3075", "1.3525", "1.3962",
    "1.4386", "1.4799", "1.5202", "1.5596", "1.5981", "1.6358",
    "1.6728", "1.7091", "1.7447", "1.7797", "1.8141", "1.8480",
    "1.8814", "1.9143", "1.9467", "1.9787", "2.0103", "2.0415",
    "2.0723", "2.1027", "2.1328", "2.1625", "2.1919", "2.2210",
    "2.2498", "2.2783", "2.3066", "2.3346", "2.3624", "2.3899",
    "2.4143",
)

COVER: Tuple[CoverInterval, ...] = tuple(
    CoverInterval(COVER_ENDPOINTS[i], COVER_ENDPOINTS[i + 1])
    for i in range(len(COVER_ENDPOINTS) - 1)
)


# ---------------------------------------------------------------------------
# Derive the six L_j polynomials with SymPy
# ---------------------------------------------------------------------------


def build_L_expressions(verbose: bool = False) -> List[sp.Expr]:
    """Return the six multi-affine polynomials L_j(alpha,C,Y,rr).

    The definitions are those of the explosion-side barrier calculation.  SymPy
    builds the expression B(t), reduces it modulo the algebraic identities for
    C,Y,rr, divides by t^2(1-t), converts the resulting quintic from the power
    basis to the Bernstein basis, and multiplies by the D_j factors.
    """

    a = alpha
    gamma = 1 + 2 * a

    # Convert scaled variables back to the variables used in the phase-plane
    # vector field formulas.
    cr = C / (1 + 3 * a)
    V2 = Y / (1 + 3 * a)
    Q2 = (cr - V2) / a
    cb = cr - 1 / (1 + 3 * a)

    # Pq and Pv are the vector-field components on the explosion side.
    Pq = Q * (
        -(V - cr) ** 2 * (V * (1 + 3 * a) - 1)
        + a * (V - cr) * (V * (V - 1) + a * Q**2)
        + a**2 * cb / gamma * Q**2
    )
    Pv = (V - cr) * (
        -V * (V - 1) * (V - cr)
        + a * Q**2 * (cr - cb / gamma - 1 + 3 * a * V)
    )

    # Quadratic barrier path F(t).  B is the barrier comparison expression.
    F = V2 * t**2 + Q2 * rr * (t - t**2)
    B = Pv.subs({V: F, Q: Q2 * t}) - sp.diff(F, t) / Q2 * Pq.subs(
        {V: F, Q: Q2 * t}
    )

    # We start working with the numerator of B.
    B_num, B_den = sp.fraction(sp.cancel(B))

    # We collect the identities solved by rr, C, Y.
    relations = [
        4 * C**2 + (6 * a - 10) * C - 60 * a**2 - 33 * a + 6,
        (4 * a + 2) * Y**2
        + (-6 * a * C - C + 6 * a**2 - a - 4) * Y
        - 2 * C**2
        + 6 * a * C
        + 5 * C,
        25 * (1 + 2 * a) * rr**2 - 6 * a,
    ]

    # We will use G to reduce factors of rr^2, Y^2 and C^2 using the relations
    # above
    G = sp.groebner(relations, rr, Y, C, order="lex", domain=sp.QQ.frac_field(a))

    # Reduce each t coefficient of B modulo the Groebner basis (one at a time).
    B_poly = sp.Poly(B_num, t)
    reduced_B_num = 0
    for k in range(B_poly.degree() + 1):
        coeff = B_poly.coeff_monomial(t**k)
        if coeff:
            reduced_B_num += G.reduce(coeff)[1] * t**k

    # B(t) is divisible by t^2(1-t).  The remaining
    # quotient Q_reduced is the degree-five polynomial whose Bernstein
    # coefficients are checked.
    Q_reduced = sp.cancel(reduced_B_num / (t**2 * (1 - t) * B_den))
    Q_num, Q_den = sp.fraction(Q_reduced)
    Q_poly = sp.Poly(Q_num, t)
    q_coeffs = [sp.cancel(Q_poly.coeff_monomial(t**i) / Q_den) for i in range(6)]

    # Compute Q(t) degree-5 Bernstein coefficients.
    b_coeffs = []
    for j in range(6):
        bj = sum(
            q_coeffs[i] * sp.binomial(j, i) / sp.binomial(5, i)
            for i in range(j + 1)
        )
        b_coeffs.append(sp.cancel(bj))

    # Positive denominator-clearing factors.
    denominators = [
        200 * a * (1 + 2 * a) ** 2 * (1 + 3 * a) ** 5,
        20000 * a**2 * (1 + 2 * a) ** 4 * (1 + 3 * a) ** 5,
        200000 * a**3 * (1 + 2 * a) ** 7 * (1 + 3 * a) ** 5,
        400000 * a**4 * (1 + 2 * a) ** 4 * (1 + 3 * a) ** 5,
        40000 * a**3 * (1 + 2 * a) ** 6 * (1 + 3 * a) ** 5,
        800 * a**2 * (1 + 2 * a) ** 6 * (1 + 3 * a) ** 5,
    ]

    Ls: List[sp.Expr] = []
    for j in range(6):
        candidate = sp.cancel(-b_coeffs[j] * denominators[j])
        num, den = sp.fraction(candidate)
        reduced_num = G.reduce(num)[1]
        L = sp.cancel(reduced_num / den)
        if sp.denom(L) != 1:
            raise ArithmeticError(f"L_{j} still has a denominator: {sp.denom(L)}")
        Ls.append(sp.expand(L))

    if verbose:
        #  This reports the size of the expression.
        for j, L in enumerate(Ls):
            poly = sp.Poly(L, alpha, C, Y, rr)
            max_alpha_degree = max(term[0][0] for term in poly.terms())
            print(
                f"normal form {j}: {len(poly.terms())} monomial terms, "
                f"max alpha degree {max_alpha_degree}"
            )

    return Ls


# ---------------------------------------------------------------------------
# Bernstein conversion
# ---------------------------------------------------------------------------

# Given a polynomial and a degree, it computes the Bernstein coefficients of
# such degree associated to the polynomial; returned as a  list of rational
# numbers
def bernstein_coeffs_from_power(
    power_coeffs: Sequence[sp.Rational], degree: Optional[int] = None
) -> List[sp.Rational]:
    """Power basis to Bernstein basis on [0,1].

    If p(u)=sum_i a_i u^i and p(u)=sum_k binom(n,k)c_k u^k(1-u)^(n-k), then

        c_k = sum_{i=0}^k a_i binom(k,i)/binom(n,i).

    """

    # The list is in power-basis order: [a_0, a_1, ...] for a_0+a_1*u+...
    polynomial_degree = len(power_coeffs) - 1
    n = polynomial_degree if degree is None else degree
    if n < polynomial_degree:
        raise ValueError("Bernstein degree cannot be smaller than polynomial degree")

    coeffs = list(power_coeffs) + [sp.Rational(0)] * (n + 1 - len(power_coeffs))
    return [
        sum(coeffs[i] * sp.binomial(k, i) / sp.binomial(n, i) for i in range(k + 1))
        for k in range(n + 1)
    ]


def restrict_to_interval_power_coeffs(
    power_coeffs: Sequence[sp.Rational], left: sp.Rational, right: sp.Rational
) -> List[sp.Rational]:
    """Return power coefficients of p(left + (right-left)u)."""

    # Substitute alpha = left + (right-left)u; so u lives on [0,1] and we can
    # check the sign of the polynomial by checking the sign of the Bernstein
    # coefficients in the variable u.
    width = right - left
    n = len(power_coeffs) - 1
    left_powers = [sp.Rational(1)] * (n + 1)
    width_powers = [sp.Rational(1)] * (n + 1)
    for i in range(1, n + 1):
        left_powers[i] = left_powers[i - 1] * left
        width_powers[i] = width_powers[i - 1] * width

    out = [sp.Rational(0)] * (n + 1)
    for m, coeff in enumerate(power_coeffs):
        if coeff == 0:
            continue
        for i in range(m + 1):
            out[i] += coeff * sp.binomial(m, i) * left_powers[m - i] * width_powers[i]

    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out

# Return the smallest Bernstein coefficient associated to a polynomial on a given
# interval [left, right]
def bernstein_min_on_interval(
    power_coeffs: Sequence[sp.Rational],
    left: sp.Rational,
    right: sp.Rational,
    degree: Optional[int] = None,
) -> sp.Rational:
    """Smallest Bernstein coefficient after restricting to [left,right]."""

    restricted = restrict_to_interval_power_coeffs(power_coeffs, left, right)
    return min(bernstein_coeffs_from_power(restricted, degree))


def bernstein_coeffs_sympy(poly: sp.Expr, var: sp.Symbol, left, right) -> List[sp.Rational]:
    """Convenience version used only for the small-alpha x-polynomials."""

    u = sp.symbols("u")
    q = sp.Poly(sp.expand(poly.subs(var, left + (right - left) * u)), u)
    n = q.degree()
    power = [q.coeff_monomial(u**i) for i in range(n + 1)]
    return bernstein_coeffs_from_power(power, n)


# ---------------------------------------------------------------------------
# Exact coefficient bookkeeping for fast cover computations
# ---------------------------------------------------------------------------


def dense_L_dictionary(L: sp.Expr) -> Dict[Tuple[int, int, int], List[sp.Rational]]:
    """Write a multi-affine L as polynomial coefficients in alpha.

    Returns a dictionary keyed by (epsilon_C, epsilon_Y, epsilon_rr).  The value
    is a dense list of coefficients in increasing powers of alpha.
    """

    # Because L is separately affine, each key only needs exponents 0 or 1
    # for C,Y,rr.  The remaining alpha-dependence is stored as a dense list.
    P = sp.Poly(sp.expand(L), alpha, C, Y, rr, domain=sp.QQ)
    sparse: Dict[Tuple[int, int, int], Dict[int, sp.Rational]] = {}
    for monomial, coeff in P.terms():
        ea, eC, eY, err = monomial
        if eC > 1 or eY > 1 or err > 1:
            raise ArithmeticError(f"L is not separately affine; saw monomial {monomial}")
        key = (eC, eY, err)
        sparse.setdefault(key, {})[ea] = coeff

    dense: Dict[Tuple[int, int, int], List[sp.Rational]] = {}
    for key, coeffs in sparse.items():
        degree = max(coeffs)
        row = [sp.Rational(0)] * (degree + 1)
        for exponent, coeff in coeffs.items():
            row[exponent] = coeff
        dense[key] = row
    return dense


def evaluate_L_power_coeffs(
    Ldict: Dict[Tuple[int, int, int], List[sp.Rational]],
    Cval: sp.Rational,
    Yval: sp.Rational,
    rrval: sp.Rational,
) -> List[sp.Rational]:
    """Substitute a box corner for C,Y,rr and return alpha-power coefficients."""

    # At a fixed box corner the C,Y,rr factors are rational constants, so L_j
    # becomes an ordinary one-variable polynomial in alpha.
    max_degree = max(len(row) - 1 for row in Ldict.values())
    out = [sp.Rational(0)] * (max_degree + 1)
    for (eC, eY, err), row in Ldict.items():
        scale = (Cval if eC else 1) * (Yval if eY else 1) * (rrval if err else 1)
        for i, coeff in enumerate(row):
            out[i] += coeff * scale

    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


# ---------------------------------------------------------------------------
# Six-decimal rational enclosures
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def C_physical(a: sp.Rational) -> sp.Expr:
    """Actual expression of C=(1+3 alpha)c_r."""

    return (sp.Rational(5) - 3 * a + sp.sqrt(1 + 102 * a + 249 * a**2)) / 4


@lru_cache(maxsize=None)
def Y_physical(a: sp.Rational) -> sp.Expr:
    """Actual expression of Y=(1+3 alpha)V_2."""

    # Y depends on the selected physical C branch.  The nested square root is
    # left symbolic until we take certified decimal floors/ceilings.
    Cval = C_physical(a)
    cr = Cval / (1 + 3 * a)
    gamma = 1 + 2 * a
    V1 = (sp.Rational(3) / (1 + 3 * a) + 2 - 2 * cr) / (3 * gamma)
    V2 = sp.Rational(1, 4) * (
        -1 + 3 * cr + 3 * V1 - sp.sqrt((1 - 3 * cr - 3 * V1) ** 2 - 24 * cr * V1)
    )
    return (1 + 3 * a) * V2


@lru_cache(maxsize=None)
def rr_magnitude_physical(a: sp.Rational) -> sp.Expr:
    """Positive magnitude of the signed paper variable rr."""

    return sp.sqrt(6 * a / (25 * (1 + 2 * a)))


def floor_decimal_expr(expr: sp.Expr, digits: int = 6) -> sp.Rational:
    """Round an exact SymPy expression downward to a rational decimal.

    Returns floor(expr * 10**digits) / 10**digits as an exact rational.
    The helper does not use numerical approximation; it verifies the floor
    by exact inequalities."""

    scale = 10**digits
    return sp.Rational(verified_scaled_floor(expr, scale), scale)


def ceil_decimal_expr(expr: sp.Expr, digits: int = 6) -> sp.Rational:
    """Round an exact SymPy expression upward to a rational decimal.

    Returns ceiling(expr * 10**digits) / 10**digits as an exact rational.
    Again, the helper does not use numerical approximation; it verifies the ceiling
    by exact inequalities."""

    scale = 10**digits
    return sp.Rational(verified_scaled_ceiling(expr, scale), scale)


def floor_decimal(q: sp.Rational, places: int = 4) -> sp.Rational:
    """Round a rational downward to a fixed number of decimal places."""

    scale = 10**places
    return sp.Rational(int(sp.floor(q * scale)), scale)


def decimal_string(q: sp.Rational, places: int = 4) -> str:
    """Decimal string rounded downward to a fixed number of places."""

    q = floor_decimal(q, places)
    whole = q.p // q.q
    rem = q.p % q.q
    digits = str((rem * (10**places)) // q.q).zfill(places)
    return f"{whole}.{digits}"


def display_decimal(q: sp.Rational) -> str:
    """Printing helper; use one decimal for very large entries, otherwise four decimals."""

    if q >= 100:
        return decimal_string(q, 1)
    return decimal_string(q, 4)


def right_endpoint_reaches_true_endpoint(interval: CoverInterval) -> bool:
    """Return True if an interval's right endpoint is at or past 1+sqrt(2)."""

    comparison = (interval.a1 - TRUE_ENDPOINT).is_nonnegative
    if comparison is None:
        raise ArithmeticError(
            f"could not certify endpoint comparison for {interval.tex_interval()}"
        )
    return bool(comparison)


def verified_scaled_floor(expr: sp.Expr, scale: int) -> int:
    """Return floor(scale*expr), verified by exact sign checks."""

    scaled = sp.simplify(scale * expr)
    candidate = sp.floor(scaled)
    if candidate.is_Integer is not True:
        raise ArithmeticError(f"could not compute exact floor of {scaled}")
    if (scaled - candidate).is_nonnegative is not True:
        raise ArithmeticError(f"uncertified floor lower inequality for {scaled}")
    if (scaled - (candidate + 1)).is_negative is not True:
        raise ArithmeticError(f"uncertified floor upper inequality for {scaled}")
    return int(candidate)


def verified_scaled_ceiling(expr: sp.Expr, scale: int) -> int:
    """Return ceiling(scale*expr), verified by exact sign checks."""

    scaled = sp.simplify(scale * expr)
    candidate = sp.ceiling(scaled)
    if candidate.is_Integer is not True:
        raise ArithmeticError(f"could not compute exact ceiling of {scaled}")
    if (scaled - candidate).is_nonpositive is not True:
        raise ArithmeticError(f"uncertified ceiling upper inequality for {scaled}")
    if (scaled - (candidate - 1)).is_positive is not True:
        raise ArithmeticError(f"uncertified ceiling lower inequality for {scaled}")
    return int(candidate)


def rational_box(
    interval: CoverInterval, digits: int = 6
) -> Tuple[sp.Rational, sp.Rational, sp.Rational, sp.Rational, sp.Rational, sp.Rational]:
    """Return a rational box containing C,Y,rr over one alpha interval; taking
    advantage of their monotonicity properties.

    C is increasing, Y is decreasing, and rr is negative and decreasing.  The
    final interval extends beyond that irrational endpoint, so its lower Y-bound
    is set to 0.
    """

    lo, hi = interval.a0, interval.a1

    # Monotonicity lets us bound C, Y, rr on the entire interval by using the
    # endpoint values. We use floors for lower bounds and ceilings for upper
    # bounds so the box is safely outward-rounded.

    Cminus = floor_decimal_expr(C_physical(lo), digits)
    Cplus = ceil_decimal_expr(C_physical(hi), digits)
    Yminus = sp.Rational(0) if right_endpoint_reaches_true_endpoint(interval) else floor_decimal_expr(Y_physical(hi), digits)
    Yplus = ceil_decimal_expr(Y_physical(lo), digits)
    rrminus = -ceil_decimal_expr(rr_magnitude_physical(hi), digits)
    rrplus = -floor_decimal_expr(rr_magnitude_physical(lo), digits)

    return Cminus, Cplus, Yminus, Yplus, rrminus, rrplus


# ---------------------------------------------------------------------------
# Small-alpha and rational-cover verification
# ---------------------------------------------------------------------------


def small_alpha_table(Ls: Sequence[sp.Expr]) -> List[Tuple[int, int, sp.Rational]]:
    """Check 0 < alpha <= 10^(-2) using x=sqrt(alpha)."""

    # The small-alpha argument is written in x=sqrt(alpha).  Factoring the first x-power leaves a
    # polynomial that can be checked by Bernstein coefficients.
    x = sp.symbols("x")
    rows: List[Tuple[int, int, sp.Rational]] = []
    for j, L in enumerate(Ls):
        minimum: Optional[sp.Rational] = None
        exponents = set()
        for Ccorner in [sp.Rational(3, 2), sp.Rational(3, 2) + 13 * x**2]:
            for Ycorner in [sp.Rational(3, 4) - 7 * x**2, sp.Rational(3, 4)]:
                for rrcorner in [-sp.Rational(1, 2) * x, -sp.Rational(12, 25) * x]:
                    # Substitute one corner of the small-alpha box.  If the
                    # resulting polynomial is positive at every such corner,
                    # we obtain positivity throughout the box.
                    expr = sp.expand(L.subs({alpha: x**2, C: Ccorner, Y: Ycorner, rr: rrcorner}))
                    p = sp.Poly(expr, x)
                    powers = [mon[0][0] for mon in p.terms() if mon[1] != 0]
                    first_power = min(powers)
                    exponents.add(first_power)
                    reduced = sp.expand(expr / x**first_power)
                    local_min = min(bernstein_coeffs_sympy(reduced, x, 0, sp.Rational(1, 10)))
                    minimum = local_min if minimum is None else min(minimum, local_min)

        if len(exponents) != 1 or minimum is None or minimum <= 0:
            raise AssertionError(f"small-alpha check failed for j={j}")
        rows.append((j, sorted(exponents)[0], minimum))
    return rows


def verify_interval(
    idx: int,
    interval: CoverInterval,
    Ldicts: Sequence[Dict[Tuple[int, int, int], List[sp.Rational]]],
    box_digits: int = 6,
    bernstein_degree: Optional[int] = None,
) -> CoverResult:
    """Verify all six L_j on one rational cover interval."""

    Cminus, Cplus, Yminus, Yplus, rrminus, rrplus = rational_box(interval, box_digits)
    minimum: Optional[sp.Rational] = None
    coefficient_j: Optional[int] = None

    # Each L_j is affine in C,Y,rr, so checking all eight corners of the box is
    # enough to bound the whole box.
    for j, Ldict in enumerate(Ldicts):
        for Ccorner, Ycorner, rrcorner in product(
            [Cminus, Cplus], [Yminus, Yplus], [rrminus, rrplus]
        ):
            power_coeffs = evaluate_L_power_coeffs(Ldict, Ccorner, Ycorner, rrcorner)
            local_min = bernstein_min_on_interval(
                power_coeffs, interval.a0, interval.a1, bernstein_degree
            )
            if minimum is None or local_min < minimum:
                minimum = local_min
                coefficient_j = j

    assert minimum is not None and coefficient_j is not None
    if minimum <= 0:
        raise AssertionError(
            f"failed on interval {interval.tex_interval()}, j={coefficient_j}, minimum={minimum}"
        )

    return CoverResult(
        idx=idx,
        interval=interval,
        minimum=minimum,
        coefficient_j=coefficient_j,
        C_bounds=(Cminus, Cplus),
        Y_bounds=(Yminus, Yplus),
        rr_bounds=(rrminus, rrplus),
    )


def verify_cover(
    Ldicts: Sequence[Dict[Tuple[int, int, int], List[sp.Rational]]],
    intervals: Sequence[CoverInterval],
    box_digits: int = 6,
    bernstein_degree: Optional[int] = None,
) -> List[CoverResult]:
    """Run the verifier over a supplied list of rational cover intervals."""

    # The caller supplies intervals.
    return [
        verify_interval(i, interval, Ldicts, box_digits, bernstein_degree)
        for i, interval in enumerate(intervals, 1)
    ]


# ---------------------------------------------------------------------------
# Reporting helpers and command-line entry point
# ---------------------------------------------------------------------------


def write_tex_tables(
    path: Path,
    small_rows: Sequence[Tuple[int, int, sp.Rational]],
    cover_rows: Sequence[CoverResult],
) -> None:
    """Write the TeX tables.

    The first table records the small-alpha exponent m_j and Bernstein lower
    bound.  The second table records the interval lower bound mu_I for each
    rational cover interval.
    """

    with path.open("w", encoding="utf-8") as handle:
        handle.write("% Generated by verify_bernstein_q_rational_cover.py\n")
        handle.write("\\begin{tabular}{c|c|c}\n")
        handle.write("$j$ & $m_j$ & smallest Bernstein coefficient of $\\mathsf S_j$ \\\\ \\hline\n")
        for j, exponent, value in small_rows:
            handle.write(f"{j} & {exponent} & ${sp.latex(value)}$ \\\\ \n")
        handle.write("\\end{tabular}\n\n")

        handle.write("\\begin{longtable}{c|c|c}\n")
        handle.write("$\\#$ & $I$ & $\\mu_I$ \\\\ \\hline\n")
        for row in cover_rows:
            handle.write(
                f"{row.idx} & $[{row.interval.left},{row.interval.right}]$ & "
                f"${display_decimal(row.minimum)}$\\\\\n"
            )
        handle.write("\\end{longtable}\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Parse CLI options, derive L_j, and run both verification passes."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--box-digits",
        type=int,
        default=6,
        help="decimal digits used in rational C/Y/rr endpoint boxes",
    )
    parser.add_argument(
        "--bernstein-degree",
        type=int,
        default=None,
        help=(
            "optional common Bernstein degree for interval polynomials; "
            "by default each polynomial uses its own degree"
        ),
    )
    parser.add_argument(
        "--tex-tables",
        type=Path,
        help="create or overwrite this file with the generated TeX tables",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print normal-form details and every cover interval",
    )
    args = parser.parse_args(argv)
    start = perf_counter()

    # Sanity-check that the cover really reaches the intended endpoints.
    if COVER[-1].a1 != RATIONAL_ENDPOINT:
        raise AssertionError("cover does not end at the declared rational endpoint")
    if not right_endpoint_reaches_true_endpoint(COVER[-1]):
        raise AssertionError("cover does not reach 1+sqrt(2)")

    # Build the symbolic L_j forms with SymPy, then convert them to list of
    # coefficients so the interval checks are performed by doing exact rational arithmetic.
    if args.verbose:
        print("Building the six L_j polynomials from the phase-plane definitions...")
    Ls = build_L_expressions(verbose=args.verbose)
    Ldicts = [dense_L_dictionary(L) for L in Ls]

    # Verify that the Groebner reduction gave the expected multi-affine shape.

    if args.verbose:
        print("\nMulti-affine check for L_j:")
    for j, L in enumerate(Ls):
        num = sp.Poly(L, alpha, C, Y, rr, domain=sp.QQ)
        degC = sp.Poly(L, C).degree()
        degY = sp.Poly(L, Y).degree()
        deg_rr = sp.Poly(L, rr).degree()
        if max(degC, degY, deg_rr) > 1:
            raise AssertionError(f"L_{j} is not affine in C,Y,rr")
        if args.verbose:
            print(f"  L_{j}: terms={len(num.terms()):3d}, deg_C={degC}, deg_Y={degY}, deg_rr={deg_rr}")

    # Run the small-alpha proof and the rational-cover proof.
    small_rows = small_alpha_table(Ls)
    cover_rows = verify_cover(Ldicts, COVER, args.box_digits, args.bernstein_degree)

    print("Small-alpha check passed.")
    print(f"Rational-cover check passed on {len(cover_rows)} intervals in [10^(-2),1+sqrt(2)].")
    print(f"Box endpoints use {args.box_digits} decimal digits; interval checks use exact rational arithmetic.")
    print(f"Global minimum lower bound in the cover: {display_decimal(min(row.minimum for row in cover_rows))}")

    if args.verbose:
        print("\nSmall-alpha rows:")
        print("j  m_j  smallest Bernstein coefficient of S_j on [0,1/10]")
        for j, exponent, value in small_rows:
            print(f"{j}  {exponent:>3}  {value}")

        print("\nUnified rational cover:")
        print("#   interval              lower bound    coefficient")
        for row in cover_rows:
            print(
                f"{row.idx:2d} [{row.interval.left},{row.interval.right}] "
                f"mu={display_decimal(row.minimum):>12} coefficient j={row.coefficient_j}"
            )

    if args.tex_tables:
        write_tex_tables(args.tex_tables, small_rows, cover_rows)
        print(f"Wrote TeX tables to {args.tex_tables}")

    print(f"Total verification time: {perf_counter() - start:.2f} seconds.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
