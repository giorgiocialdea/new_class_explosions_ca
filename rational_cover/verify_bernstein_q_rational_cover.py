#!/usr/bin/env python3
"""Exact rational-cover verification for Lemma Bernstein:q.

The script proves the sign of the six explosion-side Bernstein coefficients
B_j from Lemma~Bernstein:q.  It starts from the phase-plane definitions of
P_v, P_q and the quadratic barrier, derives the separately affine forms
in

    scaled_c = (1+3*alpha)c_r,   scaled_v = (1+3*alpha)V_2,   rr < 0,

and then checks a finite rational cover using exact rational arithmetic.
SymPy is used only to take exact floors/ceilings of the algebraic endpoint
values that define the rational boxes for (scaled_c, scaled_v, rr).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# One-variable rational-polynomial arithmetic in alpha
# ---------------------------------------------------------------------------

Poly = Tuple[Fraction, ...]  # coefficients in increasing powers of alpha
ZERO_POLY: Poly = ()
ONE_POLY: Poly = (Fraction(1),)
ALPHA_POLY: Poly = (Fraction(0), Fraction(1))
SIGMA_POLY: Poly = (Fraction(1), Fraction(3))  # 1+3 alpha
GAMMA_POLY: Poly = (Fraction(1), Fraction(2))  # 1+2 alpha


def trim(poly: Iterable[Fraction]) -> Poly:
    out = tuple(Fraction(c) for c in poly)
    n = len(out)
    while n and out[n - 1] == 0:
        n -= 1
    return out[:n]


def poly_add(p: Poly, q: Poly) -> Poly:
    n = max(len(p), len(q))
    out = [Fraction(0)] * n
    for i in range(n):
        out[i] = (p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
    return trim(out)


def poly_neg(p: Poly) -> Poly:
    return tuple(-c for c in p)


def poly_sub(p: Poly, q: Poly) -> Poly:
    return poly_add(p, poly_neg(q))


def poly_mul(p: Poly, q: Poly) -> Poly:
    if not p or not q:
        return ZERO_POLY
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, ci in enumerate(p):
        if ci == 0:
            continue
        for j, cj in enumerate(q):
            if cj:
                out[i + j] += ci * cj
    return trim(out)


def poly_scale(p: Poly, c: Fraction | int) -> Poly:
    c = Fraction(c)
    if c == 0 or not p:
        return ZERO_POLY
    return trim(ci * c for ci in p)


@lru_cache(maxsize=None)
def factor_power(name: str, exponent: int) -> Poly:
    factor = {"a": ALPHA_POLY, "s": SIGMA_POLY, "g": GAMMA_POLY}[name]
    out = ONE_POLY
    for _ in range(exponent):
        out = poly_mul(out, factor)
    return out


def poly_div_alpha_exact(p: Poly) -> Optional[Poly]:
    p = trim(p)
    if not p:
        return ZERO_POLY
    if p[0] != 0:
        return None
    return trim(p[1:])


def poly_div_linear_exact(p: Poly, factor: Poly) -> Optional[Poly]:
    """Divide by b0+b1*alpha.  Return None if the division is not exact."""

    p = trim(p)
    if not p:
        return ZERO_POLY
    if len(p) == 1:
        return None
    b0, b1 = factor
    degree = len(p) - 1
    quotient = [Fraction(0)] * degree
    remainder = list(p)
    for k in range(degree - 1, -1, -1):
        coeff = remainder[k + 1] / b1
        quotient[k] = coeff
        remainder[k + 1] -= b1 * coeff
        remainder[k] -= b0 * coeff
    if all(entry == 0 for entry in remainder):
        return trim(quotient)
    return None


@dataclass(frozen=True)
class RatPoly:
    """A rational function with denominator alpha^da(1+3a)^ds(1+2a)^dg."""

    num: Poly = ZERO_POLY
    da: int = 0
    ds: int = 0
    dg: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "num", trim(self.num))
        if not self.num:
            object.__setattr__(self, "da", 0)
            object.__setattr__(self, "ds", 0)
            object.__setattr__(self, "dg", 0)
            return
        self._cancel_denominator_factors()

    def _cancel_denominator_factors(self) -> None:
        num, da, ds, dg = self.num, self.da, self.ds, self.dg
        changed = True
        while changed and num:
            changed = False
            if da > 0:
                q = poly_div_alpha_exact(num)
                if q is not None:
                    num, da, changed = q, da - 1, True
                    continue
            if ds > 0:
                q = poly_div_linear_exact(num, SIGMA_POLY)
                if q is not None:
                    num, ds, changed = q, ds - 1, True
                    continue
            if dg > 0:
                q = poly_div_linear_exact(num, GAMMA_POLY)
                if q is not None:
                    num, dg, changed = q, dg - 1, True
                    continue
        object.__setattr__(self, "num", num)
        object.__setattr__(self, "da", da)
        object.__setattr__(self, "ds", ds)
        object.__setattr__(self, "dg", dg)

    @staticmethod
    def const(value: Fraction | int) -> "RatPoly":
        value = Fraction(value)
        return RatPoly((value,) if value else ZERO_POLY)

    @staticmethod
    def alpha() -> "RatPoly":
        return RatPoly(ALPHA_POLY)

    def is_zero(self) -> bool:
        return not self.num

    def __neg__(self) -> "RatPoly":
        return RatPoly(poly_neg(self.num), self.da, self.ds, self.dg)

    def __add__(self, other: "RatPoly | Fraction | int") -> "RatPoly":
        other = to_ratpoly(other)
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        da = max(self.da, other.da)
        ds = max(self.ds, other.ds)
        dg = max(self.dg, other.dg)
        n1, n2 = self.num, other.num
        if da > self.da:
            n1 = poly_mul(n1, factor_power("a", da - self.da))
        if ds > self.ds:
            n1 = poly_mul(n1, factor_power("s", ds - self.ds))
        if dg > self.dg:
            n1 = poly_mul(n1, factor_power("g", dg - self.dg))
        if da > other.da:
            n2 = poly_mul(n2, factor_power("a", da - other.da))
        if ds > other.ds:
            n2 = poly_mul(n2, factor_power("s", ds - other.ds))
        if dg > other.dg:
            n2 = poly_mul(n2, factor_power("g", dg - other.dg))
        return RatPoly(poly_add(n1, n2), da, ds, dg)

    __radd__ = __add__

    def __sub__(self, other: "RatPoly | Fraction | int") -> "RatPoly":
        return self + (-to_ratpoly(other))

    def __rsub__(self, other: "RatPoly | Fraction | int") -> "RatPoly":
        return to_ratpoly(other) + (-self)

    def __mul__(self, other: "RatPoly | Fraction | int") -> "RatPoly":
        other = to_ratpoly(other)
        if self.is_zero() or other.is_zero():
            return RatPoly()
        return RatPoly(poly_mul(self.num, other.num), self.da + other.da, self.ds + other.ds, self.dg + other.dg)

    __rmul__ = __mul__

    def scale(self, value: Fraction | int) -> "RatPoly":
        return RatPoly(poly_scale(self.num, Fraction(value)), self.da, self.ds, self.dg)

    def div_alpha(self, exponent: int = 1) -> "RatPoly":
        return RatPoly(self.num, self.da + exponent, self.ds, self.dg)

    def div_sigma(self, exponent: int = 1) -> "RatPoly":
        return RatPoly(self.num, self.da, self.ds + exponent, self.dg)

    def div_gamma(self, exponent: int = 1) -> "RatPoly":
        return RatPoly(self.num, self.da, self.ds, self.dg + exponent)

    def mul_alpha(self, exponent: int = 1) -> "RatPoly":
        return RatPoly(poly_mul(self.num, factor_power("a", exponent)), self.da, self.ds, self.dg)

    def mul_sigma(self, exponent: int = 1) -> "RatPoly":
        return RatPoly(poly_mul(self.num, factor_power("s", exponent)), self.da, self.ds, self.dg)

    def mul_gamma(self, exponent: int = 1) -> "RatPoly":
        return RatPoly(poly_mul(self.num, factor_power("g", exponent)), self.da, self.ds, self.dg)

    def as_poly_after_multiplication(self, const: int, alpha_power: int, sigma_power: int, gamma_power: int) -> Poly:
        out = RatPoly(poly_scale(self.num, Fraction(const)), self.da, self.ds, self.dg)
        if alpha_power:
            out = out.mul_alpha(alpha_power)
        if sigma_power:
            out = out.mul_sigma(sigma_power)
        if gamma_power:
            out = out.mul_gamma(gamma_power)
        if out.da or out.ds or out.dg:
            raise ArithmeticError(f"remaining denominator alpha^{out.da} sigma^{out.ds} gamma^{out.dg}")
        return out.num


def to_ratpoly(value: RatPoly | Fraction | int) -> RatPoly:
    if isinstance(value, RatPoly):
        return value
    return RatPoly.const(value)


# ---------------------------------------------------------------------------
# Normal-form arithmetic in scaled_c, scaled_v, and rr<0
# ---------------------------------------------------------------------------

BasisKey = Tuple[int, int, int]  # powers of scaled_c, scaled_v, rr
NormalFormDict = Dict[BasisKey, RatPoly]
PolynomialDict = Dict[BasisKey, Poly]

# scaled_c^2 = ((60a^2+33a-6) + (10-6a) scaled_c)/4.
SCALED_C_SQUARE: NormalFormDict = {
    (0, 0, 0): RatPoly((Fraction(-3, 2), Fraction(33, 4), Fraction(15))),
    (1, 0, 0): RatPoly((Fraction(5, 2), Fraction(-3, 2))),
}

# rr^2 = 6a/(25(1+2a)).
RR_SQUARE = RatPoly(ALPHA_POLY, dg=1).scale(Fraction(6, 25))

# The scaled_v^2 relation, with scaled_c^2 already reduced.
SCALED_V_SQUARE: NormalFormDict = {
    (0, 1, 0): RatPoly((Fraction(4), Fraction(1), Fraction(-6)), dg=1).scale(Fraction(1, 2)),
    (1, 1, 0): RatPoly((Fraction(1), Fraction(6)), dg=1).scale(Fraction(1, 2)),
    (0, 0, 0): SCALED_C_SQUARE[(0, 0, 0)].scale(2).div_gamma().scale(Fraction(1, 2)),
    (1, 0, 0): (SCALED_C_SQUARE[(1, 0, 0)].scale(2) - RatPoly((Fraction(5), Fraction(6)))).div_gamma().scale(Fraction(1, 2)),
}


@lru_cache(maxsize=None)
def reduce_monomial(c_power: int, v_power: int, rr_power: int) -> NormalFormDict:
    if rr_power >= 2:
        return {key: value * RR_SQUARE for key, value in reduce_monomial(c_power, v_power, rr_power - 2).items()}
    if c_power >= 2:
        out: NormalFormDict = {}
        for (dc, dv, dr), coeff in SCALED_C_SQUARE.items():
            for key, value in reduce_monomial(c_power - 2 + dc, v_power + dv, rr_power + dr).items():
                out[key] = out.get(key, RatPoly()) + coeff * value
        return {key: value for key, value in out.items() if not value.is_zero()}
    if v_power >= 2:
        out = {}
        for (dc, dv, dr), coeff in SCALED_V_SQUARE.items():
            for key, value in reduce_monomial(c_power + dc, v_power - 2 + dv, rr_power + dr).items():
                out[key] = out.get(key, RatPoly()) + coeff * value
        return {key: value for key, value in out.items() if not value.is_zero()}
    return {(c_power, v_power, rr_power): RatPoly.const(1)}


class NormalForm:
    __slots__ = ("data",)

    def __init__(self, data: Optional[NormalFormDict] = None):
        self.data: NormalFormDict = {}
        if data:
            for key, value in data.items():
                value = to_ratpoly(value)
                if not value.is_zero():
                    self.data[key] = value

    @staticmethod
    def const(value: RatPoly | Fraction | int) -> "NormalForm":
        return NormalForm({(0, 0, 0): to_ratpoly(value)})

    def is_zero(self) -> bool:
        return not self.data

    def __neg__(self) -> "NormalForm":
        return NormalForm({key: -value for key, value in self.data.items()})

    def __add__(self, other: "NormalForm | RatPoly | Fraction | int") -> "NormalForm":
        other = to_normalform(other)
        out = dict(self.data)
        for key, value in other.data.items():
            out[key] = out.get(key, RatPoly()) + value
        return NormalForm(out)

    __radd__ = __add__

    def __sub__(self, other: "NormalForm | RatPoly | Fraction | int") -> "NormalForm":
        return self + (-to_normalform(other))

    def __rsub__(self, other: "NormalForm | RatPoly | Fraction | int") -> "NormalForm":
        return to_normalform(other) + (-self)

    def __mul__(self, other: "NormalForm | RatPoly | Fraction | int") -> "NormalForm":
        other = to_normalform(other)
        out: NormalFormDict = {}
        for (c1, v1, r1), coeff1 in self.data.items():
            for (c2, v2, r2), coeff2 in other.data.items():
                for key, reduction_coeff in reduce_monomial(c1 + c2, v1 + v2, r1 + r2).items():
                    out[key] = out.get(key, RatPoly()) + coeff1 * coeff2 * reduction_coeff
        return NormalForm(out)

    __rmul__ = __mul__

    def scale(self, value: Fraction | int) -> "NormalForm":
        return NormalForm({key: coeff.scale(value) for key, coeff in self.data.items()})

    def div_alpha(self) -> "NormalForm":
        return NormalForm({key: coeff.div_alpha() for key, coeff in self.data.items()})

    def div_sigma(self) -> "NormalForm":
        return NormalForm({key: coeff.div_sigma() for key, coeff in self.data.items()})

    def div_gamma(self) -> "NormalForm":
        return NormalForm({key: coeff.div_gamma() for key, coeff in self.data.items()})

    def as_polynomial_dict_after_multiplication(self, const: int, alpha_power: int, sigma_power: int, gamma_power: int) -> PolynomialDict:
        return {
            key: coeff.as_poly_after_multiplication(const, alpha_power, sigma_power, gamma_power)
            for key, coeff in self.data.items()
        }


def to_normalform(value: NormalForm | RatPoly | Fraction | int) -> NormalForm:
    if isinstance(value, NormalForm):
        return value
    return NormalForm.const(value)


ZERO_NF = NormalForm()
ONE_NF = NormalForm.const(1)
ALPHA_NF = NormalForm.const(RatPoly.alpha())
SCALED_C_NF = NormalForm({(1, 0, 0): RatPoly.const(1)})
SCALED_V_NF = NormalForm({(0, 1, 0): RatPoly.const(1)})
RR_NF = NormalForm({(0, 0, 1): RatPoly.const(1)})
SIGMA_NF = NormalForm.const(RatPoly(SIGMA_POLY))


# ---------------------------------------------------------------------------
# Polynomial arithmetic in t with normal-form coefficients
# ---------------------------------------------------------------------------

NFPoly = List[NormalForm]


def nfpoly_trim(poly: NFPoly) -> NFPoly:
    while poly and poly[-1].is_zero():
        poly.pop()
    return poly


def nfpoly_add(p: NFPoly, q: NFPoly) -> NFPoly:
    n = max(len(p), len(q))
    out = []
    for i in range(n):
        out.append((p[i] if i < len(p) else ZERO_NF) + (q[i] if i < len(q) else ZERO_NF))
    return nfpoly_trim(out)


def nfpoly_neg(p: NFPoly) -> NFPoly:
    return [-entry for entry in p]


def nfpoly_sub(p: NFPoly, q: NFPoly) -> NFPoly:
    return nfpoly_add(p, nfpoly_neg(q))


def nfpoly_mul(p: NFPoly, q: NFPoly) -> NFPoly:
    if not p or not q:
        return []
    out = [ZERO_NF for _ in range(len(p) + len(q) - 1)]
    for i, pi in enumerate(p):
        if pi.is_zero():
            continue
        for j, qj in enumerate(q):
            if qj.is_zero():
                continue
            out[i + j] = out[i + j] + pi * qj
    return nfpoly_trim(out)


def nfpoly_shift(p: NFPoly, shift: int) -> NFPoly:
    return [ZERO_NF] * shift + p


def nfpoly_scale(p: NFPoly, scalar: NormalForm | RatPoly | Fraction | int) -> NFPoly:
    scalar = to_normalform(scalar)
    return [entry * scalar for entry in p]


def nfpoly_derivative(p: NFPoly) -> NFPoly:
    return [p[i].scale(i) for i in range(1, len(p))]


def nfpoly_power(p: NFPoly, exponent: int) -> NFPoly:
    out = [ONE_NF]
    for _ in range(exponent):
        out = nfpoly_mul(out, p)
    return out


def build_normal_forms(verbose: bool = False) -> List[PolynomialDict]:
    """Derive the six normal forms from the phase-plane definitions."""

    c_r = SCALED_C_NF.div_sigma()
    V_2 = SCALED_V_NF.div_sigma()
    Q_2 = (SCALED_C_NF - SCALED_V_NF).div_alpha().div_sigma()
    c_b = (SCALED_C_NF - ONE_NF).div_sigma()

    # F(t)=V_2 t^2 + Q_2 rr(t-t^2), with rr=rr<0.
    F_poly = nfpoly_add(
        nfpoly_shift([V_2], 2),
        nfpoly_mul(nfpoly_shift([Q_2 * RR_NF], 1), [ONE_NF, -ONE_NF]),
    )
    Q_poly = nfpoly_shift([Q_2], 1)

    V_minus_cr = nfpoly_sub(F_poly, [c_r])
    sigma_V_minus_one = nfpoly_sub(nfpoly_scale(F_poly, SIGMA_NF), [ONE_NF])
    V_minus_one = nfpoly_sub(F_poly, [ONE_NF])
    Q_square = nfpoly_power(Q_poly, 2)

    bracket_q = nfpoly_add(
        nfpoly_add(
            nfpoly_neg(nfpoly_mul(nfpoly_power(V_minus_cr, 2), sigma_V_minus_one)),
            nfpoly_scale(
                nfpoly_mul(V_minus_cr, nfpoly_add(nfpoly_mul(F_poly, V_minus_one), nfpoly_scale(Q_square, ALPHA_NF))),
                ALPHA_NF,
            ),
        ),
        nfpoly_scale(Q_square, NormalForm.const(RatPoly(ALPHA_POLY).mul_alpha().div_gamma()) * c_b),
    )

    # Since P_q=Q_2 t * bracket_q, the quotient P_q/Q_2 is t*bracket_q.
    Pq_over_Q2 = nfpoly_shift(bracket_q, 1)

    Pv_inner = nfpoly_add(
        nfpoly_neg(nfpoly_mul(nfpoly_mul(F_poly, V_minus_one), V_minus_cr)),
        nfpoly_scale(Q_square, ALPHA_NF * (c_r - c_b.div_gamma() - ONE_NF)),
    )
    Pv_inner = nfpoly_add(
        Pv_inner,
        nfpoly_scale(nfpoly_mul(Q_square, F_poly), NormalForm.const(RatPoly(ALPHA_POLY).mul_alpha().scale(3))),
    )
    Pv = nfpoly_mul(V_minus_cr, Pv_inner)
    B_poly = nfpoly_sub(Pv, nfpoly_mul(nfpoly_derivative(F_poly), Pq_over_Q2))

    # B(t)=t^2(1-t)Q(t).  Recover the six power coefficients of Q recursively.
    q_coefficients: List[NormalForm] = []
    for j in range(6):
        coeff = B_poly[j + 2] if j + 2 < len(B_poly) else ZERO_NF
        if j:
            coeff = coeff + q_coefficients[-1]
        q_coefficients.append(coeff)

    if len(B_poly) > 8 and any(not entry.is_zero() for entry in B_poly[9:]):
        raise ArithmeticError("unexpected degree after division by t^2(1-t)")
    high_check = (B_poly[8] if len(B_poly) > 8 else ZERO_NF) + q_coefficients[-1]
    if not high_check.is_zero():
        raise ArithmeticError("division by t^2(1-t) failed in the top coefficient")

    # d_j(alpha) values are const*alpha^p*(1+3alpha)^q*(1+2alpha)^r.
    denominators = [
        (200, 1, 5, 2),
        (20000, 2, 5, 4),
        (200000, 3, 5, 7),
        (400000, 4, 5, 4),
        (40000, 3, 5, 6),
        (800, 2, 5, 6),
    ]

    normal_forms: List[PolynomialDict] = []
    for j in range(6):
        bernstein_coeff = ZERO_NF
        for i in range(j + 1):
            bernstein_coeff = bernstein_coeff + q_coefficients[i].scale(Fraction(comb(j, i), comb(5, i)))
        const, alpha_power, sigma_power, gamma_power = denominators[j]
        normal_form = (-bernstein_coeff).as_polynomial_dict_after_multiplication(
            const, alpha_power, sigma_power, gamma_power
        )
        normal_forms.append(normal_form)

    if verbose:
        for j, form in enumerate(normal_forms):
            max_degree = max(len(poly) for poly in form.values()) - 1
            print(f"normal form {j}: {len(form)} monomial terms, max alpha degree {max_degree}")
    return normal_forms


# ---------------------------------------------------------------------------
# Bernstein conversion and rational-cover checks
# ---------------------------------------------------------------------------

alpha_symbol = sp.symbols("alpha")
TRUE_ENDPOINT = 1 + sp.sqrt(2)
RATIONAL_ENDPOINT = Fraction(24143, 10000)


@dataclass(frozen=True)
class CoverInterval:
    left: str
    right: str

    @property
    def a0(self) -> Fraction:
        return Fraction(self.left)

    @property
    def a1(self) -> Fraction:
        return Fraction(self.right)

    def tex_interval(self) -> str:
        return f"[{self.left},{self.right}]"


BASE_COVER: Tuple[CoverInterval, ...] = tuple(
    CoverInterval(a, b)
    for a, b in [
        ("0.0100", "0.0476"),
        ("0.0476", "0.0997"),
        ("0.0997", "0.1233"),
        ("0.1233", "0.1436"),
        ("0.1436", "0.1638"),
        ("0.1638", "0.1858"),
        ("0.1858", "0.2116"),
        ("0.2116", "0.2435"),
        ("0.2435", "0.2765"),
        ("0.2765", "0.3047"),
        ("0.3047", "0.3294"),
        ("0.3294", "0.3517"),
        ("0.3517", "0.3726"),
        ("0.3726", "0.3919"),
        ("0.3919", "0.4091"),
        ("0.4091", "0.4247"),
        ("0.4247", "0.4394"),
        ("0.4394", "0.4537"),
        ("0.4537", "0.4679"),
        ("0.4679", "0.4823"),
        ("0.4823", "0.4971"),
        ("0.4971", "0.5124"),
        ("0.5124", "0.5285"),
        ("0.5285", "0.5455"),
        ("0.5455", "0.5546"),
        ("0.5546", "0.5637"),
        ("0.5637", "0.5832"),
        ("0.5832", "0.6042"),
        ("0.6042", "0.6269"),
        ("0.6269", "0.6516"),
        ("0.6516", "0.6785"),
        ("0.6785", "0.7078"),
        ("0.7078", "0.7398"),
        ("0.7398", "0.7748"),
        ("0.7748", "0.8131"),
        ("0.8131", "0.8549"),
        ("0.8549", "0.9006"),
        ("0.9006", "0.9504"),
        ("0.9504", "1.0000"),
    ]
)

EXTENSION_ENDPOINTS = [
    "1.0000", "1.0569", "1.1111", "1.1630", "1.2129", "1.2610", "1.3075",
    "1.3525", "1.3962", "1.4386", "1.4799", "1.5202", "1.5596", "1.5981",
    "1.6358", "1.6728", "1.7091", "1.7447", "1.7797", "1.8141", "1.8480",
    "1.8814", "1.9143", "1.9467", "1.9787", "2.0103", "2.0415", "2.0723",
    "2.1027", "2.1328", "2.1625", "2.1919", "2.2210", "2.2498", "2.2783",
    "2.3066", "2.3346", "2.3624", "2.3899", "2.4143",
]

EXTENSION_COVER: Tuple[CoverInterval, ...] = tuple(
    CoverInterval(EXTENSION_ENDPOINTS[i], EXTENSION_ENDPOINTS[i + 1])
    for i in range(len(EXTENSION_ENDPOINTS) - 1)
)


def full_cover() -> Tuple[CoverInterval, ...]:
    return BASE_COVER + EXTENSION_COVER


def fraction_to_sympy(q: Fraction) -> sp.Rational:
    return sp.Rational(q.numerator, q.denominator)


def fraction_from_sympy(q: sp.Rational) -> Fraction:
    return Fraction(int(q.p), int(q.q))


@lru_cache(maxsize=None)
def scaled_c_value(a: Fraction) -> sp.Expr:
    x = fraction_to_sympy(a)
    return (sp.Rational(5) - 3 * x + sp.sqrt(1 + 102 * x + 249 * x**2)) / 4


@lru_cache(maxsize=None)
def scaled_v_value(a: Fraction) -> sp.Expr:
    x = fraction_to_sympy(a)
    scaled_c = scaled_c_value(a)
    c_r = scaled_c / (1 + 3 * x)
    gamma = 1 + 2 * x
    V_1 = (sp.Rational(3) / (1 + 3 * x) + 2 - 2 * c_r) / (3 * gamma)
    V_2 = sp.Rational(1, 4) * (
        -1 + 3 * c_r + 3 * V_1 - sp.sqrt((1 - 3 * c_r - 3 * V_1) ** 2 - 24 * c_r * V_1)
    )
    return (1 + 3 * x) * V_2


@lru_cache(maxsize=None)
def abs_rr_value(a: Fraction) -> sp.Expr:
    x = fraction_to_sympy(a)
    return sp.sqrt(6 * x / (25 * (1 + 2 * x)))


def exact_floor_decimal(expr: sp.Expr, digits: int) -> Fraction:
    scale = 10**digits
    value = sp.floor(expr * scale)
    if not value.is_Integer:
        raise ArithmeticError(f"could not certify floor for {expr}")
    return Fraction(int(value), scale)


def exact_ceil_decimal(expr: sp.Expr, digits: int) -> Fraction:
    scale = 10**digits
    value = sp.ceiling(expr * scale)
    if not value.is_Integer:
        raise ArithmeticError(f"could not certify ceiling for {expr}")
    return Fraction(int(value), scale)


def interval_exceeds_true_endpoint(interval: CoverInterval) -> bool:
    comparison = (fraction_to_sympy(interval.a1) - TRUE_ENDPOINT).is_positive
    if comparison is None:
        return bool(sp.N(fraction_to_sympy(interval.a1) - TRUE_ENDPOINT, 80) > 0)
    return bool(comparison)


@lru_cache(maxsize=None)
def scaled_c_floor(a: Fraction, digits: int) -> Fraction:
    return exact_floor_decimal(scaled_c_value(a), digits)


@lru_cache(maxsize=None)
def scaled_c_ceil(a: Fraction, digits: int) -> Fraction:
    return exact_ceil_decimal(scaled_c_value(a), digits)


@lru_cache(maxsize=None)
def scaled_v_floor(a: Fraction, digits: int) -> Fraction:
    return exact_floor_decimal(scaled_v_value(a), digits)


@lru_cache(maxsize=None)
def scaled_v_ceil(a: Fraction, digits: int) -> Fraction:
    return exact_ceil_decimal(scaled_v_value(a), digits)


@lru_cache(maxsize=None)
def abs_rr_floor(a: Fraction, digits: int) -> Fraction:
    return exact_floor_decimal(abs_rr_value(a), digits)


@lru_cache(maxsize=None)
def abs_rr_ceil(a: Fraction, digits: int) -> Fraction:
    return exact_ceil_decimal(abs_rr_value(a), digits)


def rational_box(interval: CoverInterval, digits: int) -> Tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
    a0, a1 = interval.a0, interval.a1
    scaled_c_min = scaled_c_floor(a0, digits)
    scaled_c_max = scaled_c_ceil(a1, digits)
    scaled_v_max = scaled_v_ceil(a0, digits)
    if interval_exceeds_true_endpoint(interval):
        scaled_v_min = Fraction(0)
    else:
        scaled_v_min = scaled_v_floor(a1, digits)
    # rr is negative and decreases with alpha.
    rr_min = -abs_rr_ceil(a1, digits)
    rr_max = -abs_rr_floor(a0, digits)
    return scaled_c_min, scaled_c_max, scaled_v_min, scaled_v_max, rr_min, rr_max


def evaluate_normal_form_power_coefficients(form: PolynomialDict, scaled_c: Fraction, scaled_v: Fraction, rr: Fraction) -> List[Fraction]:
    max_degree = max(len(poly) for poly in form.values())
    out = [Fraction(0)] * max_degree
    for (c_power, v_power, rr_power), poly in form.items():
        scale = (scaled_c if c_power else 1) * (scaled_v if v_power else 1) * (rr if rr_power else 1)
        for i, coeff in enumerate(poly):
            out[i] += coeff * scale
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def restrict_power_coefficients(poly: Sequence[Fraction], left: Fraction, right: Fraction) -> List[Fraction]:
    width = right - left
    n = len(poly) - 1
    left_powers = [Fraction(1)] * (n + 1)
    width_powers = [Fraction(1)] * (n + 1)
    for i in range(1, n + 1):
        left_powers[i] = left_powers[i - 1] * left
        width_powers[i] = width_powers[i - 1] * width
    out = [Fraction(0)] * (n + 1)
    for m, coeff in enumerate(poly):
        if coeff == 0:
            continue
        for i in range(m + 1):
            out[i] += coeff * Fraction(comb(m, i)) * left_powers[m - i] * width_powers[i]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def bernstein_coefficients_from_power(poly: Sequence[Fraction], degree: Optional[int] = None) -> List[Fraction]:
    m = len(poly) - 1
    n = m if degree is None else degree
    if n < m:
        raise ValueError("Bernstein degree cannot be smaller than polynomial degree")
    coeffs = list(poly) + [Fraction(0)] * (n + 1 - len(poly))
    return [
        sum(coeffs[i] * Fraction(comb(k, i), comb(n, i)) for i in range(k + 1))
        for k in range(n + 1)
    ]


def bernstein_minimum_on_interval(poly: Sequence[Fraction], left: Fraction, right: Fraction, degree: Optional[int] = None) -> Fraction:
    restricted = restrict_power_coefficients(poly, left, right)
    return min(bernstein_coefficients_from_power(restricted, degree))


@dataclass(frozen=True)
class CoverResult:
    idx: int
    interval: CoverInterval
    minimum: Fraction
    coefficient_j: int
    box: Tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]


def verify_interval(
    idx: int,
    interval: CoverInterval,
    normal_forms: Sequence[PolynomialDict],
    box_digits: int,
    bernstein_degree: Optional[int] = None,
) -> CoverResult:
    scaled_c_min, scaled_c_max, scaled_v_min, scaled_v_max, rr_min, rr_max = rational_box(interval, box_digits)
    minimum: Optional[Fraction] = None
    coefficient_j: Optional[int] = None
    for j, form in enumerate(normal_forms):
        for c_corner, v_corner, rr_corner in product(
            [scaled_c_min, scaled_c_max], [scaled_v_min, scaled_v_max], [rr_min, rr_max]
        ):
            power_coefficients = evaluate_normal_form_power_coefficients(form, c_corner, v_corner, rr_corner)
            local_minimum = bernstein_minimum_on_interval(
                power_coefficients, interval.a0, interval.a1, bernstein_degree
            )
            if minimum is None or local_minimum < minimum:
                minimum = local_minimum
                coefficient_j = j
    assert minimum is not None and coefficient_j is not None
    if minimum <= 0:
        raise AssertionError(
            f"failed on interval {interval.tex_interval()}, coefficient j={coefficient_j}, minimum={minimum}"
        )
    return CoverResult(
        idx=idx,
        interval=interval,
        minimum=minimum,
        coefficient_j=coefficient_j,
        box=(scaled_c_min, scaled_c_max, scaled_v_min, scaled_v_max, rr_min, rr_max),
    )


def decimal_down(value: Fraction, digits: int = 4) -> str:
    scale = 10**digits
    scaled = value.numerator * scale // value.denominator
    sign = "-" if scaled < 0 else ""
    scaled = abs(scaled)
    return f"{sign}{scaled // scale}.{scaled % scale:0{digits}d}"


def display_decimal(value: Fraction) -> str:
    return decimal_down(value, 1 if value >= 100 else 4)


# ---------------------------------------------------------------------------
# Small-alpha verification
# ---------------------------------------------------------------------------


def poly_x_add(p: List[Fraction], q: List[Fraction]) -> List[Fraction]:
    n = max(len(p), len(q))
    out = [Fraction(0)] * n
    for i in range(n):
        out[i] = (p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_x_mul(p: List[Fraction], q: List[Fraction]) -> List[Fraction]:
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, ci in enumerate(p):
        if ci == 0:
            continue
        for j, cj in enumerate(q):
            if cj:
                out[i + j] += ci * cj
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def alpha_poly_to_x_squared(poly: Poly) -> List[Fraction]:
    out = [Fraction(0)] * (2 * (len(poly) - 1) + 1)
    for i, coeff in enumerate(poly):
        out[2 * i] = coeff
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def small_alpha_table(normal_forms: Sequence[PolynomialDict]) -> List[Tuple[int, int, Fraction]]:
    c_corners = [[Fraction(3, 2)], [Fraction(3, 2), Fraction(0), Fraction(13)]]
    v_corners = [[Fraction(3, 4), Fraction(0), Fraction(-7)], [Fraction(3, 4)]]
    rr_corners = [[Fraction(0), Fraction(-1, 2)], [Fraction(0), Fraction(-12, 25)]]
    rows: List[Tuple[int, int, Fraction]] = []
    for j, form in enumerate(normal_forms):
        minimum: Optional[Fraction] = None
        exponents = set()
        for c_corner in c_corners:
            for v_corner in v_corners:
                for rr_corner in rr_corners:
                    total = [Fraction(0)]
                    for (c_power, v_power, rr_power), alpha_poly in form.items():
                        term = alpha_poly_to_x_squared(alpha_poly)
                        if c_power:
                            term = poly_x_mul(term, c_corner)
                        if v_power:
                            term = poly_x_mul(term, v_corner)
                        if rr_power:
                            term = poly_x_mul(term, rr_corner)
                        total = poly_x_add(total, term)
                    first_nonzero = next(i for i, coeff in enumerate(total) if coeff)
                    exponents.add(first_nonzero)
                    reduced = total[first_nonzero:]
                    local_minimum = bernstein_minimum_on_interval(reduced, Fraction(0), Fraction(1, 10))
                    if minimum is None or local_minimum < minimum:
                        minimum = local_minimum
        if len(exponents) != 1 or minimum is None or minimum <= 0:
            raise AssertionError(f"small-alpha check failed for j={j}")
        rows.append((j, next(iter(exponents)), minimum))
    return rows


def verify_cover(
    normal_forms: Sequence[PolynomialDict],
    intervals: Sequence[CoverInterval],
    box_digits: int,
    bernstein_degree: Optional[int] = None,
) -> List[CoverResult]:
    return [
        verify_interval(i, interval, normal_forms, box_digits, bernstein_degree)
        for i, interval in enumerate(intervals, 1)
    ]


def write_tex_tables(path: Path, small_rows: Sequence[Tuple[int, int, Fraction]], cover_rows: Sequence[CoverResult]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("% Generated by verify_bernstein_q_rational_cover.py\n")
        handle.write("\\begin{tabular}{c|c|c}\n")
        handle.write("$j$ & $m_j$ & smallest Bernstein coefficient of $\\mathsf S_j$ \\\\ \\hline\n")
        for j, exponent, value in small_rows:
            handle.write(f"{j} & {exponent} & ${sp.latex(fraction_to_sympy(value))}$ \\\\ \n")
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--box-digits", type=int, default=6, help="decimal digits used in rational endpoint boxes")
    parser.add_argument(
        "--bernstein-degree",
        type=int,
        default=None,
        help="optional common Bernstein degree for interval polynomials; by default each polynomial uses its own degree",
    )
    parser.add_argument(
        "--tex-tables",
        type=Path,
        help="create or overwrite this file with the generated TeX tables",
    )
    parser.add_argument("--verbose", action="store_true", help="print normal-form details and every interval")
    args = parser.parse_args(argv)

    normal_forms = build_normal_forms(verbose=args.verbose)
    small_rows = small_alpha_table(normal_forms)
    cover_rows = verify_cover(normal_forms, full_cover(), args.box_digits, args.bernstein_degree)

    print("Small-alpha check passed.")
    print(f"Rational-cover check passed on {len(cover_rows)} intervals in [10^(-2),1+sqrt(2)].")
    print(f"Box endpoints use {args.box_digits} decimal digits and exact rational arithmetic.")
    print(f"Global minimum lower bound in the cover: {display_decimal(min(row.minimum for row in cover_rows))}")

    if args.verbose:
        for row in cover_rows:
            print(
                f"{row.idx:2d} [{row.interval.left},{row.interval.right}] "
                f"mu={display_decimal(row.minimum):>12} coefficient j={row.coefficient_j}"
            )

    if args.tex_tables:
        write_tex_tables(args.tex_tables, small_rows, cover_rows)
        print(f"Wrote TeX tables to {args.tex_tables}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
