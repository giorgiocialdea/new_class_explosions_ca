
r"""Lemma-level interval verification for Table~\ref{table:bb}.

The routine ``bernstein_bounds()`` corresponds to the paper statement
Lemma~\ref{lemma:Bernstein:q} / Table~\ref{table:bb}.  The file is
self-contained apart from the parameters in ``parameters.py``
and the constants in ``supplementary_data/bounds.txt``.
"""

from __future__ import annotations

import time
from pathlib import Path
from dataclasses import dataclass, field
from decimal import Decimal
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from sage.all import QQ, RealBallField, RealIntervalField

import parameters
from parameters import (
    B_INTERVAL,
    CUTS,
    MAX_BOXES,
    MAX_DEPTH,
    PRECISION_BITS,
    PROGRESS_EVERY,
)

########################################
#Set the verbosity from parameters.py

def set_verbose(level: int = 1, *, print_subintervals=None) -> None:
    """Set runtime verbosity from script.py."""
    parameters.VERBOSE = int(level)
    if print_subintervals is not None:
        parameters.PRINT_SUBINTERVALS = bool(print_subintervals)

########################################
#Sanity check to verify no RBF has unresolved uncertainty.

def o(x: Any) -> bool:
    """Return True when Sage did not print an uncertainty marker for x."""
    return "?" not in str(x)

########################################
# Progress printer (if the verbosity requires it)
def print_iter_N(index, total=None, text="") -> None:
    if not parameters.VERBOSE:
        return
    if total is None:
        print(f"  [{index}] {text}", flush=True)
    else:
        print(f"  [{index}/{total}] {text}", flush=True)

########################################
# Import the bounds from supplementary_data/bounds.txt.
def load_bounds(path=None):
    """Load supplementary_data/bounds.txt.

    Each non-comment line has the form

        label|value|side|description

    The numerical value is converted to an RBF element, and the original string
    and metadata are also returned, indexed by label.
    """
    if path is None:
        path = Path(__file__).resolve().parent / "supplementary_data" / "bounds.txt"
    else:
        path = Path(path)

    bounds = {}
    bounds_str = {}
    metadata = {}

    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [part.strip() for part in line.split("|")]
            if len(parts) < 4:
                raise ValueError(
                    f"bad bounds.txt line {line_no}: expected "
                    "label|value|side|description"
                )
            label, value, side = parts[:3]
            description = "|".join(parts[3:]).strip()
            if side not in {"lower", "upper"}:
                raise ValueError(
                    f"bad bounds.txt line {line_no}: side must be 'lower' or 'upper'"
                )
            if label in bounds:
                raise ValueError(f"duplicate bounds.txt label {label!r}")
            bounds[label] = parameters.RBF(value)
            bounds_str[label] = value
            metadata[label] = {
                "line": line_no,
                "side": side,
                "description": description,
            }

    return bounds, bounds_str, metadata

#####################################
# Detailed printer for accepted boxes.  It reports the local interval
# enclosure and the difference from the requested lower/upper table bounds.
def _print_local_box(problem, box, value, method, lower_bound, upper_bound) -> None:
    """Verbose debugging output for one accepted subinterval.

    """
    lo = lower(value)
    hi = upper(value)
    lower_threshold = _bound_lower_threshold(lower_bound)
    upper_threshold = _bound_upper_threshold(upper_bound)
    lower_margin = lo - lower_threshold
    upper_margin = upper_threshold - hi
    lower_rel = lower_margin / (1 + abs(lower_threshold))
    upper_rel = upper_margin / (1 + abs(upper_threshold))
    print(
        "    "
        f"{problem.name} I=[{qq_to_str(box.a)}, {qq_to_str(box.b)}] "
        f"depth={box.depth} method={method}; "
        f"enclosure=[{lo}, {hi}]; "
        f"margins=(lower {lower_margin}, upper {upper_margin}); "
        f"relative_margins=(lower {lower_rel}, upper {upper_rel})",
        flush=True,
    )

#################################################
# Load all table constants.  BOUNDS contains Sage RBF values;
# BOUNDS_STR keeps the original decimal strings for printing output
# BOUNDS_META contains the metadata relative to the bound.
BOUNDS, BOUNDS_STR, BOUNDS_META = load_bounds()

################################################
# To avoid accidental rounding from Python floats
# all endpoints and constants are converted exactly to rationals.
def to_QQ(x: Any) -> Any:
    """Convert int / Decimal / 'a/b' string / decimal string to Sage QQ exactly.
    """
    if hasattr(x, "parent") and x.parent() == QQ:
        return x
    if isinstance(x, int):
        return QQ(x)
    if isinstance(x, Decimal):
        n, den = x.as_integer_ratio()
        return QQ(n) / QQ(den)
    if isinstance(x, float):
        raise TypeError(
            "Python float input is forbidden. Use a string like '0.1' or a fraction '1/10'."
        )
    s = str(x).strip()
    if "/" in s:
        return QQ(s)
    d = Decimal(s)
    n, den = d.as_integer_ratio()
    return QQ(n) / QQ(den)

#################################################
# Convert a Sage rational back to text for the printing output.
def qq_to_str(q: Any) -> str:
    return str(QQ(q))


# ---------------------------------------------------------------------------
# Ball arithmetic context (with per-evaluation cache for Shared nodes)
# ---------------------------------------------------------------------------
#################################
# Context that owns the RBF parameters and a cache for shared subexpressions.
class Context:
    """Holds the working RealBallField/RealIntervalField and a per-pass cache.

    The cache is cleared at the start of every top-level eval_direct or
    eval_derivative call.
    """

    def __init__(self, prec: int):
        self.prec = int(prec)
        self.RBF = RealBallField(self.prec)
        self.RIF = RealIntervalField(self.prec)
        self.cache: Dict[Any, Any] = {}

    def reset_cache(self) -> None:
        self.cache.clear()

    def ball_interval(self, a: Any, b: Any) -> Any:
        a = to_QQ(a)
        b = to_QQ(b)
        if b < a:
            raise ValueError(f"bad interval [{a},{b}]")
        return self.RBF(self.RIF(a, b))

    def ball_point(self, x: Any) -> Any:
        return self.RBF(to_QQ(x))

###################################
# Convenience wrappers that return Sage interval endpoints.
def lower(x: Any) -> Any:
    return x.lower()


def upper(x: Any) -> Any:
    return x.upper()

###################################
# Helpers that check the sign of a given RBF.
def strictly_positive(x: Any) -> bool:
    return bool(lower(x) > 0)


def nonnegative(x: Any) -> bool:
    return bool(lower(x) >= 0)




##############################################
# In this class formulas are converted into expression trees; each expression
# is evaluated using interval arithmetic and can be automatically
# differentiated.
class Expr:
    def __add__(self, other): return Add(self, as_expr(other))
    def __radd__(self, other): return Add(as_expr(other), self)
    def __sub__(self, other): return Add(self, -as_expr(other))
    def __rsub__(self, other): return Add(as_expr(other), -self)
    def __mul__(self, other): return Mul(self, as_expr(other))
    def __rmul__(self, other): return Mul(as_expr(other), self)
    def __truediv__(self, other): return Div(self, as_expr(other))
    def __rtruediv__(self, other): return Div(as_expr(other), self)
    def __neg__(self): return Mul(Const(QQ(-1)), self)

    def __pow__(self, n: int) -> "Expr":
        if not isinstance(n, int):
            raise TypeError("Only nonnegative integer powers are supported.")
        if n < 0:
            raise ValueError("Only nonnegative integer powers are supported.")
        return Pow(self, n)

    def eval_ad(self, x, ctx, want_derivative):
        raise NotImplementedError

##############################################
# Constant node in the expression tree.
#Transforms a given rational number into a node of the expression tree,
#treating the constant as an exact RBF element. If want_derivative, it returns
# zero (that is the derivative of a constant is zero)
@dataclass(frozen=True)
class Const(Expr):
    q: Any

    def __post_init__(self):
        object.__setattr__(self, "q", to_QQ(self.q))

    def eval_ad(self, x, ctx, want_derivative):
        val = ctx.RBF(self.q)
        der = ctx.RBF(0) if want_derivative else None
        return True, val, der, ""

    def __str__(self):
        return qq_to_str(self.q)

##############################################
# Variable node in the expression tree. Given the RBF element x
# as an input it returns x. If want_derivative, it returns 1.
@dataclass(frozen=True)
class Var(Expr):
    name: str

    def eval_ad(self, x, ctx, want_derivative):
        der = ctx.RBF(1) if want_derivative else None
        return True, x, der, ""

    def __str__(self):
        return self.name

##############################################
# Addition node in the expression tree.
# The value and derivatives are combined term-by-term.
@dataclass(frozen=True)
class Add(Expr):
    left: Expr
    right: Expr

    def eval_ad(self, x, ctx, want_derivative):
        ok1, v1, d1, r1 = self.left.eval_ad(x, ctx, want_derivative)
        if not ok1:
            return False, None, None, r1
        ok2, v2, d2, r2 = self.right.eval_ad(x, ctx, want_derivative)
        if not ok2:
            return False, None, None, r2
        val = v1 + v2
        der = d1 + d2 if want_derivative else None
        return True, val, der, ""

    def __str__(self):
        return f"({self.left} + {self.right})"

##############################################
# Multiplication node in the expression tree.
# The derivative is computed from Leibniz rule.
@dataclass(frozen=True)
class Mul(Expr):
    left: Expr
    right: Expr

    def eval_ad(self, x, ctx, want_derivative):
        ok1, v1, d1, r1 = self.left.eval_ad(x, ctx, want_derivative)
        if not ok1:
            return False, None, None, r1
        ok2, v2, d2, r2 = self.right.eval_ad(x, ctx, want_derivative)
        if not ok2:
            return False, None, None, r2
        val = v1 * v2
        der = d1 * v2 + v1 * d2 if want_derivative else None
        return True, val, der, ""

    def __str__(self):
        return f"({self.left}*{self.right})"

##############################################
# Division node in the expression tree.
# Before dividing the denominator is verified
# away from zero
@dataclass(frozen=True)
class Div(Expr):
    numerator: Expr
    denominator: Expr

    def eval_ad(self, x, ctx, want_derivative):
        ok1, v1, d1, r1 = self.numerator.eval_ad(x, ctx, want_derivative)
        if not ok1:
            return False, None, None, r1
        ok2, v2, d2, r2 = self.denominator.eval_ad(x, ctx, want_derivative)
        if not ok2:
            return False, None, None, r2
        if not (strictly_positive(v2) or bool(upper(v2) < 0)):
            return False, None, None, "division denominator not verified away from zero"
        val = v1 / v2
        der = (d1 * v2 - v1 * d2) / (v2 ** 2) if want_derivative else None
        return True, val, der, ""

    def __str__(self):
        return f"({self.numerator}/{self.denominator})"

##############################################
# Power node in the expression tree.
# Only allows nonnegative integer powers
@dataclass(frozen=True)
class Pow(Expr):
    base: Expr
    n: int

    def eval_ad(self, x, ctx, want_derivative):
        ok, v, d, r = self.base.eval_ad(x, ctx, want_derivative)
        if not ok:
            return False, None, None, r
        if self.n == 0:
            val = ctx.RBF(1)
            der = ctx.RBF(0) if want_derivative else None
        else:
            val = v ** self.n
            der = ctx.RBF(self.n) * (v ** (self.n - 1)) * d if want_derivative else None
        return True, val, der, ""

    def __str__(self):
        return f"({self.base}**{self.n})"

##############################################
# Square-root node. For evaluation, the radicand
# must be nonnegative. For derivative evaluation, it must be strictly positive.
@dataclass(frozen=True)
class Sqrt(Expr):
    arg: Expr
    label: Optional[str] = None

    def eval_ad(self, x, ctx, want_derivative):
        ok, rad, drad, reason = self.arg.eval_ad(x, ctx, want_derivative)
        if not ok:
            return False, None, None, reason
        lbl = self.label or "sqrt"
        if want_derivative:
            if not strictly_positive(rad):
                return False, None, None, (
                    f"sqrt radicand '{lbl}' not verified strictly positive for derivative"
                )
        else:
            if not nonnegative(rad):
                return False, None, None, f"sqrt radicand '{lbl}' not verified nonnegative"
        val = rad.sqrt()
        der = drad / (ctx.RBF(2) * val) if want_derivative else None
        return True, val, der, ""

    def __str__(self):
        if self.label:
            return f"sqrt[{self.label}]({self.arg})"
        return f"sqrt({self.arg})"

##############################################
# This node is used to cache repeated subexpressions
# like cr or V2.

class Shared(Expr):
    """Memoizing wrapper.
    Each (id(self), want_derivative) is evaluated at most once per box.
    because the cache lives on Context and is cleared at the start of every
    top-level eval_direct/eval_derivative call.
    """

    __slots__ = ("inner", "label")

    def __init__(self, inner: Expr, label: Optional[str] = None):
        self.inner = inner
        self.label = label

    def eval_ad(self, x, ctx, want_derivative):
        key = (id(self), bool(want_derivative))
        hit = ctx.cache.get(key)
        if hit is not None:
            return hit
        result = self.inner.eval_ad(x, ctx, want_derivative)
        ctx.cache[key] = result
        return result

    def __str__(self):
        return f"<{self.label}>" if self.label else str(self.inner)


ExprLike = Union[Expr, int, str]


##############################################
# Helper that makes sure constants are converted to
# the class Expr
def as_expr(x: Any) -> Expr:
    if isinstance(x, Expr):
        return x
    return Const(x)

##############################################
# Helper to construct sqrt nodes

def sqrt_expr(x: Any, label: Optional[str] = None) -> Expr:
    return Sqrt(as_expr(x), label=label)


# ---------------------------------------------------------------------------
# Problem API used by config files
# ---------------------------------------------------------------------------


########################################
#To build expressions
class ProblemAPI:
    def __init__(self, default_interval: Optional[Sequence[Any]] = None):
        # The DSL is one-variable.  `alpha` is kept as the historical name,
        # while `b` is the name used by the b5b* configuration.
        self.alpha = Var("alpha")
        self.b = Var("b")
        self.default_interval = (
            tuple(default_interval) if default_interval is not None else None
        )

    def const(self, q):
        return Const(q)


    def sqrt(self, x, label=None):
        return sqrt_expr(x, label=label)

    def shared(self, x, label=None):
        """Wrap a subexpression so its interval value is cached per box.
        """
        e = as_expr(x)
        if isinstance(e, Shared):
            return e
        return Shared(e, label=label)


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------
##########################
# Result object returned by direct, derivative, and Taylor evaluation routines.
@dataclass
class EvalResult:
    ok: bool
    value: Optional[Any] = None
    method: str = "direct"
    reason: str = ""

# Direct interval evaluation on either a point ball or an interval box.
def eval_direct(expr: Expr, x: Any, ctx: Context) -> EvalResult:
    ctx.reset_cache()
    ok, val, _, reason = expr.eval_ad(x, ctx, want_derivative=False)
    return EvalResult(ok=ok, value=val, method="direct", reason=reason)

# Evaluate a rigorous interval enclosure of the derivative on an interval.
def eval_derivative(expr: Expr, x: Any, ctx: Context) -> EvalResult:
    ctx.reset_cache()
    ok, _, der, reason = expr.eval_ad(x, ctx, want_derivative=True)
    return EvalResult(ok=ok, value=der, method="derivative", reason=reason)

# First-order Taylor enclosure: f([a,b]) is contained in f(mid) + ([-rho,rho])f'([a,b]).

def eval_taylor1(expr: Expr, a: Any, b: Any, ctx: Context) -> EvalResult:
    a = to_QQ(a)
    b = to_QQ(b)
    X = ctx.ball_interval(a, b)
    m = (a + b) / 2
    rho = (b - a) / 2

    mid = eval_direct(expr, ctx.ball_point(m), ctx)
    if not mid.ok:
        return EvalResult(False, method="taylor1",
                          reason="midpoint evaluation failed: " + mid.reason)

    der = eval_derivative(expr, X, ctx)
    if not der.ok:
        return EvalResult(False, method="taylor1",
                          reason="derivative enclosure failed: " + der.reason)

    delta = ctx.ball_interval(-rho, rho)
    val = mid.value + delta * der.value
    return EvalResult(True, value=val, method="taylor1")



# ---------------------------------------------------------------------------
# Adaptive proof of the table bounds
# ---------------------------------------------------------------------------
# One closed subinterval of the original domain.  The depth records how many bisections
# produced it.
@dataclass(frozen=True)
class Box:
    a: Any
    b: Any
    depth: int = 0

    def bisect(self) -> Tuple["Box", "Box"]:
        m = (self.a + self.b) / 2
        return Box(self.a, m, self.depth + 1), Box(m, self.b, self.depth + 1)

##############################################
# Verification problem for one table row.
# Contains one expression, the lower/upper bounds from bounds.txt,
# and an interval.
@dataclass(frozen=True)
class BoundProblem:
    name: str
    symbol: str
    expr: Expr
    interval: Tuple[Any, Any]
    cuts: Tuple[Any, ...]
    lower_label: str
    upper_label: str

##############################################
# Build the starting boxes by splitting the full interval at the configured cuts
# (the cuts are configured in parameters). The cuts are not necessary and they
# are just chosen for convenience/speed.
def initial_boxes(interval: Sequence[Any], cuts: Sequence[Any]) -> List[Box]:
    a = to_QQ(interval[0])
    b = to_QQ(interval[1])
    if not a < b:
        raise ValueError(f"bad interval [{a},{b}]")
    pts = sorted({a, b, *(to_QQ(c) for c in cuts)})
    if pts[0] != a or pts[-1] != b:
        raise ValueError("cuts must lie inside the problem interval")
    return [Box(pts[i], pts[i + 1], 0) for i in range(len(pts) - 1)]

##############################################
# Conservative table-bound interpretation.  For a lower bound, use the upper
# endpoint of its ball so the final inequality is strict. Similarly, for an
# upper bound use its lower bound.
def _bound_lower_threshold(bound):
    """Strict lower-table bound, interpreted safely from an RBF element."""
    return upper(bound)


def _bound_upper_threshold(bound):
    """Strict upper-table bound, interpreted safely from an RBF element."""
    return lower(bound)

##############################################
# Return True only when a computed enclosure lies strictly between the
# requested lower and upper bounds from bounds.txt.

def _inside_table_bounds(value, lower_bound, upper_bound) -> bool:
    return bool(
        lower(value) > _bound_lower_threshold(lower_bound)
        and upper(value) < _bound_upper_threshold(upper_bound)
    )
##############################################
# Try to prove one box using direct interval arithmetic; or
# if that fails, Taylor enclosure.
def _best_table_enclosure(
    expr: Expr,
    box: Box,
    ctx: Context,
    lower_bound,
    upper_bound,
) -> EvalResult:

    # Method 1: evaluate the expression directly on the entire interval box.

    X = ctx.ball_interval(box.a, box.b)
    direct = eval_direct(expr, X, ctx)
    if (
        direct.ok
        and direct.value is not None
        and _inside_table_bounds(direct.value, lower_bound, upper_bound)
    ):
        return direct

    # Method 2: if direct evaluation failed, try a first-order
    # Taylor enclosure: f(midpoint) + [-rho, rho] * f'([a,b]).
    taylor = eval_taylor1(expr, box.a, box.b, ctx)
    if (
        taylor.ok
        and taylor.value is not None
        and _inside_table_bounds(taylor.value, lower_bound, upper_bound)
    ):
        return taylor

    # Neither method proved the desired bounds, return the failure.
    if taylor.ok and taylor.value is not None:
        return EvalResult(
            False,
            value=taylor.value,
            method="taylor1",
            reason="table bounds not verified on this box",
        )
    if direct.ok and direct.value is not None:
        return EvalResult(
            False,
            value=direct.value,
            method="direct",
            reason="table bounds not verified on this box",
        )
    return direct

######################################################
# Main loop that checks the bounds for a single Bernstein coefficient.
def prove_table_row(
    problem: BoundProblem,
    *,
    lower_bound,
    upper_bound,
    prec: int = PRECISION_BITS,
    max_depth: int = MAX_DEPTH,
    max_boxes: int = MAX_BOXES,
    progress_every: int = PROGRESS_EVERY,
    print_subintervals: Optional[bool] = None,
) -> Dict[str, Any]:
    """Verify one row of Table~\\ref{table:bb} by adaptive subdivision."""
    if progress_every == 0 and parameters.VERBOSE:
        progress_every = parameters.VERBOSE_COUNTER_I
    if print_subintervals is None:
        print_subintervals = parameters.PRINT_SUBINTERVALS or parameters.VERBOSE >= 2
    ctx = Context(prec)
    stack = list(reversed(initial_boxes(problem.interval, problem.cuts)))
    failures: List[Dict[str, Any]] = []

    accepted_count = 0
    bisected_count = 0
    boxes_seen = 0
    max_depth_reached = 0
    global_lower = None
    global_upper = None
    method_counts: Dict[str, int] = {}
    t0 = time.perf_counter()

    #  A box is accepted if one enclosure method proves it;
    #  otherwise it is split into two smaller boxes.

    while stack:
        if accepted_count + len(stack) > max_boxes:
            failures.append({"reason": "maximum number of boxes exceeded"})
            break

        box = stack.pop()
        boxes_seen += 1
        max_depth_reached = max(max_depth_reached, box.depth)

        if progress_every and boxes_seen % progress_every == 0:
            print(
                f"  {problem.name}: seen={boxes_seen}, accepted={accepted_count}, "
                f"stack={len(stack)}, depth={max_depth_reached}",
                flush=True,
            )

        result = _best_table_enclosure(
            problem.expr, box, ctx, lower_bound, upper_bound
        )

        if result.ok and result.value is not None:
            # This box is rigorously inside the desired table bounds.
            accepted_count += 1
            method_counts[result.method] = method_counts.get(result.method, 0) + 1
            lo = lower(result.value)
            hi = upper(result.value)
            if global_lower is None or lo < global_lower:
                global_lower = lo
            if global_upper is None or hi > global_upper:
                global_upper = hi
            if print_subintervals:
                _print_local_box(problem, box, result.value, result.method, lower_bound, upper_bound)
            continue

        if box.depth >= max_depth:
            # The interval is still inconclusive, but the allowed subdivision
            # depth has been exhausted.
            failures.append({
                "a": qq_to_str(box.a),
                "b": qq_to_str(box.b),
                "depth": box.depth,
                "method": result.method,
                "reason": result.reason,
                "lower": str(lower(result.value)) if result.value is not None else None,
                "upper": str(upper(result.value)) if result.value is not None else None,
            })
            continue
        # The current box did not prove the bound yet, so bisect and retry each half
        left, right = box.bisect()
        stack.append(right)
        stack.append(left)
        bisected_count += 1

    verified = (not failures) and accepted_count > 0
    out: Dict[str, Any] = {
        "name": problem.name,
        "symbol": problem.symbol,
        "verified": verified,
        "box_count": accepted_count,
        "boxes_seen": boxes_seen,
        "bisected_count": bisected_count,
        "max_depth_reached": max_depth_reached,
        "method_counts": method_counts,
        "elapsed_seconds": time.perf_counter() - t0,
        "failures": failures,
    }
    if verified:
        out["global_lower_bound"] = global_lower
        out["global_upper_bound"] = global_upper
    return out



# ---------------------------------------------------------------------------
# Explicit formulae for the six rescaled Bernstein coefficients
# ---------------------------------------------------------------------------
# Shared quantities used by every b5b* coefficient formula. We use the node
# shared so in each subinterval the enclosures of such quantities
# are cached
def make_cr_and_V2(alpha, api):
    """Build cr(alpha) and V2(alpha).

    Both are wrapped with api.shared(...) so the verifier evaluates each
    only once per box.
    """
    a = alpha
    sqrt = api.sqrt
    shared = api.shared

    cr = shared(
        (5 - 3 * a + sqrt(1 + 102 * a + 249 * a**2, label="cr radicand"))
        / (4 + 12 * a),
        label="cr",
    )

    V2_rad = (
        16
        - 32 * cr
        + 17 * cr**2
        + 8 * a
        - 174 * cr * a
        + 146 * cr**2 * a
        - 47 * a**2
        - 330 * cr * a**2
        + 453 * cr**2 * a**2
        - 12 * a**3
        - 360 * cr * a**3
        + 612 * cr**2 * a**3
        + 36 * a**4
        - 216 * cr * a**4
        + 324 * cr**2 * a**4
    )

    V2 = shared(
        (
            4
            + cr
            + a
            + 9 * cr * a
            - 6 * a**2
            + 18 * cr * a**2
            - sqrt(V2_rad, label="V2 radicand")
        )
        / (4 * (1 + 2 * a) * (1 + 3 * a)),
        label="V2",
    )

    return cr, V2



###############################
# The coefficient \BB_0. Also alpha=b^2 and radicand
# sqrt(6 + 12 * alpha) are cached.

def b5b0(b, api):
    alpha = api.shared(b**2, label="alpha")
    cr, V2 = make_cr_and_V2(alpha, api)
    sqrt_alpha = b
    sqrt_6_12alpha = api.shared(
        api.sqrt(6 + 12 * alpha, label="6 + 12 alpha")
    )

    return (
        32 * cr**3 * sqrt_alpha * (1 + 3 * alpha)
        - V2**2 * sqrt_alpha * (119 + 132 * alpha)
        + cr**2
        * (
            5 * sqrt_6_12alpha
            + 15 * alpha * sqrt_6_12alpha
            - sqrt_alpha * (119 + 132 * alpha + 64 * V2 * (1 + 3 * alpha))
        )
        + cr
        * V2
        * (
            -5 * sqrt_6_12alpha
            - 15 * alpha * sqrt_6_12alpha
            + sqrt_alpha
            * (263 + 32 * V2 * (1 + 3 * alpha) + alpha * (389 + 150 * alpha))
        )
    )

###############################
# The coefficient \BB_1.
def b5b1(b, api):
    alpha = api.shared(b**2, label="alpha")
    cr, V2 = make_cr_and_V2(alpha, api)
    sqrt_alpha = b
    sqrt_6 = api.shared(api.sqrt(6, label="6"))
    sqrt_1_2alpha = api.shared(api.sqrt(1 + 2 * alpha, label="1 + 2 alpha"))

    return (
        2 * sqrt_6 * V2**3 * sqrt_alpha * (47 + 66 * alpha)
        + 5
        * cr**4
        * (1 + 3 * alpha)
        * (
            -11 * sqrt_6 * sqrt_alpha
            - 20 * sqrt_6 * alpha * sqrt_alpha
            + 6 * sqrt_1_2alpha
            + 246 * alpha * sqrt_1_2alpha
        )
        + cr**3
        * (
            2 * sqrt_6 * (28 + 95 * V2) * sqrt_alpha
            + sqrt_6 * (618 + 1145 * V2) * alpha * sqrt_alpha
            + 75 * sqrt_6 * (12 + 29 * V2) * alpha**2 * sqrt_alpha
            + 1350 * sqrt_6 * V2 * alpha**3 * sqrt_alpha
            - 60 * (-1 + V2) * sqrt_1_2alpha
            - 120 * (29 + 22 * V2) * alpha * sqrt_1_2alpha
            - 90 * (47 + 82 * V2) * alpha**2 * sqrt_1_2alpha
        )
        + cr**2
        * V2
        * (
            2 * sqrt_6 * (91 - 95 * V2) * sqrt_alpha
            - sqrt_6 * (179 + 1145 * V2) * alpha * sqrt_alpha
            - 75 * sqrt_6 * (13 + 29 * V2) * alpha**2 * sqrt_alpha
            - 450 * sqrt_6 * (1 + 3 * V2) * alpha**3 * sqrt_alpha
            + 30 * (-4 + V2) * sqrt_1_2alpha
            + 30 * (257 + 44 * V2) * alpha * sqrt_1_2alpha
            + 30 * (407 + 123 * V2) * alpha**2 * sqrt_1_2alpha
            + 4500 * alpha**3 * sqrt_1_2alpha
        )
        + cr
        * V2**2
        * (
            60 * sqrt_1_2alpha
            - 3480 * alpha * sqrt_1_2alpha
            - 4230 * alpha**2 * sqrt_1_2alpha
            + 5
            * sqrt_6
            * V2
            * sqrt_alpha
            * (1 + 3 * alpha)
            * (11 + 20 * alpha)
            + sqrt_6
            * sqrt_alpha
            * (-332 + alpha * (-571 + 75 * alpha * (1 + 6 * alpha)))
        )
    )

###############################
# The coefficient \BB_2.
def b5b2(b, api):
    alpha = api.shared(b**2, label="alpha")
    cr, V2 = make_cr_and_V2(alpha, api)
    sqrt_alpha = b
    sqrt_6 = api.shared(api.sqrt(6, label="6"))
    sqrt_1_2alpha = api.shared(api.sqrt(1 + 2 * alpha, label="1 + 2 alpha"))

    return -(
        2
        * cr**4
        * (1 + 3 * alpha)
        * (
            -30 * sqrt_6
            + 265 * sqrt_6 * alpha
            + 1775 * sqrt_6 * alpha**2
            + 2250 * sqrt_6 * alpha**3
            - 114 * sqrt_alpha * sqrt_1_2alpha
            - 9225 * alpha * sqrt_alpha * sqrt_1_2alpha
            - 17850 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
        )
        + V2**3
        * (
            30 * sqrt_6
            - 1945 * sqrt_6 * alpha
            - 6920 * sqrt_6 * alpha**2
            - 5820 * sqrt_6 * alpha**3
            + 6 * (-25 + 62 * V2) * sqrt_alpha * sqrt_1_2alpha
            + 3 * (-675 + 572 * V2) * alpha * sqrt_alpha * sqrt_1_2alpha
            + 1800 * (-3 + V2) * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            - 3900 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
        )
        + cr
        * V2**2
        * (
            -V2
            * (1 + 3 * alpha)
            * (
                -60 * sqrt_6
                + 530 * sqrt_6 * alpha
                + 3550 * sqrt_6 * alpha**2
                + 4500 * sqrt_6 * alpha**3
                + 1788 * sqrt_alpha * sqrt_1_2alpha
                + 6625 * alpha * sqrt_alpha * sqrt_1_2alpha
                + 8500 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
                + 2500 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            )
            - 5
            * (1 + 2 * alpha)
            * (
                18 * sqrt_6
                - 1353 * sqrt_6 * alpha
                - 2221 * sqrt_6 * alpha**2
                + 475 * sqrt_6 * alpha**3
                + 1650 * sqrt_6 * alpha**4
                + 180 * sqrt_alpha * sqrt_1_2alpha
                - 9535 * alpha * sqrt_alpha * sqrt_1_2alpha
                - 12850 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
                - 250 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
                + 1500 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            )
        )
        + cr**2
        * V2
        * (
            5
            * (1 + 2 * alpha)
            * (
                18 * sqrt_6
                - 978 * sqrt_6 * alpha
                - 346 * sqrt_6 * alpha**2
                + 2725 * sqrt_6 * alpha**3
                + 1650 * sqrt_6 * alpha**4
                + 450 * sqrt_alpha * sqrt_1_2alpha
                - 19410 * alpha * sqrt_alpha * sqrt_1_2alpha
                - 31905 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
                - 11250 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            )
            + 3
            * V2
            * (1 + 3 * alpha)
            * (
                -60 * sqrt_6
                + 655 * sqrt_6 * alpha
                + 5425 * sqrt_6 * alpha**2
                + 10500 * sqrt_6 * alpha**3
                + 5500 * sqrt_6 * alpha**4
                + 744 * sqrt_alpha * sqrt_1_2alpha
                - 2125 * alpha * sqrt_alpha * sqrt_1_2alpha
                - 4150 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
                + 7500 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
                + 5000 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            )
        )
        + cr**3
        * (
            -1200 * sqrt_alpha * sqrt_1_2alpha
            - 588 * V2 * sqrt_alpha * sqrt_1_2alpha
            + 40575 * alpha * sqrt_alpha * sqrt_1_2alpha
            + 29711 * V2 * alpha * sqrt_alpha * sqrt_1_2alpha
            + 141300 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            + 157325 * V2 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            + 110700 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            + 186200 * V2 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            - 7500 * V2 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            - 5
            * sqrt_6
            * (1 + 2 * alpha)
            * (
                6
                + alpha * (-26 + 3 * alpha * (431 + 750 * alpha))
                + 3
                * V2
                * (1 + 3 * alpha)
                * (-12 + 5 * alpha * (31 + 5 * alpha * (31 + 22 * alpha)))
            )
        )
    )


###############################
# The coefficient \BB_3.
def b5b3(b, api):
    alpha = api.shared(b**2, label="alpha")
    cr, V2 = make_cr_and_V2(alpha, api)
    sqrt_alpha = b
    sqrt_6 = api.shared(api.sqrt(6, label="6"))
    sqrt_1_2alpha = api.shared(api.sqrt(1 + 2 * alpha, label="1 + 2 alpha"))

    return (
        2
        * cr**4
        * (1 + 3 * alpha)
        * (
            60 * sqrt_6 * sqrt_alpha
            - 300 * sqrt_6 * alpha * sqrt_alpha
            - 1875 * sqrt_6 * alpha**2 * sqrt_alpha
            + 18 * sqrt_1_2alpha
            + 75 * alpha * sqrt_1_2alpha
            + 12325 * alpha**2 * sqrt_1_2alpha
        )
        + cr**2
        * V2
        * (
            90 * sqrt_6 * (-2 + 3 * V2) * sqrt_alpha
            + 5 * sqrt_6 * (1610 - 753 * V2) * alpha * sqrt_alpha
            + 10 * sqrt_6 * (377 - 3810 * V2) * alpha**2 * sqrt_alpha
            - 375 * sqrt_6 * (51 + 245 * V2) * alpha**3 * sqrt_alpha
            - 11250 * sqrt_6 * (1 + 5 * V2) * alpha**4 * sqrt_alpha
            + 216 * V2 * sqrt_1_2alpha
            - 72 * (50 + 41 * V2) * alpha * sqrt_1_2alpha
            + 50 * (2503 - 3 * V2) * alpha**2 * sqrt_1_2alpha
            + 25 * (8489 + 428 * V2) * alpha**3 * sqrt_1_2alpha
            + 3750 * (19 - 23 * V2) * alpha**4 * sqrt_1_2alpha
            - 67500 * V2 * alpha**5 * sqrt_1_2alpha
        )
        + V2**3
        * (
            -30 * sqrt_6 * (2 + V2) * sqrt_alpha
            + 45 * sqrt_6 * (80 - 17 * V2) * alpha * sqrt_alpha
            + 5 * sqrt_6 * (1143 - 605 * V2) * alpha**2 * sqrt_alpha
            - 500 * sqrt_6 * (3 + 7 * V2) * alpha**3 * sqrt_alpha
            - 1500 * sqrt_6 * (3 + V2) * alpha**4 * sqrt_alpha
            + 36 * V2 * sqrt_1_2alpha
            - 642 * V2 * alpha * sqrt_1_2alpha
            + 100 * (41 - 10 * V2) * alpha**2 * sqrt_1_2alpha
            + 25 * (267 + 350 * V2) * alpha**3 * sqrt_1_2alpha
            + 2500 * (1 + 8 * V2) * alpha**4 * sqrt_1_2alpha
            + 7500 * (1 + 2 * V2) * alpha**5 * sqrt_1_2alpha
        )
        + cr
        * V2**2
        * (
            -30 * sqrt_6 * (-6 + V2) * sqrt_alpha
            + 5 * sqrt_6 * (-2160 + 607 * V2) * alpha * sqrt_alpha
            + 5 * sqrt_6 * (-3354 + 4075 * V2) * alpha**2 * sqrt_alpha
            + 125 * sqrt_6 * (51 + 328 * V2) * alpha**3 * sqrt_alpha
            + 750 * sqrt_6 * (21 + 32 * V2) * alpha**4 * sqrt_alpha
            - 144 * V2 * sqrt_1_2alpha
            + 24 * (75 + 107 * V2) * alpha * sqrt_1_2alpha
            + 25 * (-2699 + 515 * V2) * alpha**2 * sqrt_1_2alpha
            + 375 * (-278 + 11 * V2) * alpha**3 * sqrt_1_2alpha
            - 1250 * (3 + 28 * V2) * alpha**4 * sqrt_1_2alpha
            - 7500 * (-3 + 5 * V2) * alpha**5 * sqrt_1_2alpha
        )
        + cr**3
        * (
            -144 * V2 * sqrt_1_2alpha
            + 1800 * alpha * sqrt_1_2alpha
            + 768 * V2 * alpha * sqrt_1_2alpha
            - 54275 * alpha**2 * sqrt_1_2alpha
            - 39325 * V2 * alpha**2 * sqrt_1_2alpha
            - 72150 * alpha**3 * sqrt_1_2alpha
            - 125025 * V2 * alpha**3 * sqrt_1_2alpha
            + 11250 * V2 * alpha**4 * sqrt_1_2alpha
            + 5
            * sqrt_6
            * sqrt_alpha
            * (
                12
                + alpha * (-170 + alpha * (1457 + 2850 * alpha))
                + V2
                * (1 + 3 * alpha)
                * (-66 + 5 * alpha * (109 + 25 * alpha * (29 + 18 * alpha)))
            )
        )
    )


###############################
# The coefficient \BB_4.
def b5b4(b, api):
    alpha = api.shared(b**2, label="alpha")
    cr, V2 = make_cr_and_V2(alpha, api)
    sqrt_alpha = b
    sqrt_6 = api.shared(api.sqrt(6, label="6"))
    sqrt_1_2alpha = api.shared(api.sqrt(1 + 2 * alpha, label="1 + 2 alpha"))

    return -(
        2
        * cr**4
        * (1 + 3 * alpha)
        * (
            -6 * sqrt_6
            + 5
            * alpha
            * (7 * sqrt_6 + 55 * sqrt_6 * alpha - 327 * sqrt_alpha * sqrt_1_2alpha)
        )
        + V2**3
        * (
            90 * sqrt_alpha * sqrt_1_2alpha
            + 210 * V2 * sqrt_alpha * sqrt_1_2alpha
            - 605 * alpha * sqrt_alpha * sqrt_1_2alpha
            + 545 * V2 * alpha * sqrt_alpha * sqrt_1_2alpha
            - 1500 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            - 3005 * V2 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            - 1000 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            - 11250 * V2 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            - 3000 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            - 9000 * V2 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            + sqrt_6
            * V2
            * (1 + 2 * alpha)
            * (1 + 3 * alpha)
            * (-18 + 25 * alpha * (7 + 6 * alpha))
            + sqrt_6
            * (6 + alpha * (-613 + 12 * alpha * (-89 + 50 * alpha * (1 + 3 * alpha))))
        )
        + cr**2
        * V2
        * (
            570 * sqrt_alpha * sqrt_1_2alpha
            + 330 * V2 * sqrt_alpha * sqrt_1_2alpha
            - 16095 * alpha * sqrt_alpha * sqrt_1_2alpha
            - 865 * V2 * alpha * sqrt_alpha * sqrt_1_2alpha
            - 28290 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            - 1565 * V2 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            - 9000 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            + 16500 * V2 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            + 13500 * V2 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            + sqrt_6
            * V2
            * (1 + 3 * alpha)
            * (-90 + alpha * (727 + 725 * alpha * (7 + 6 * alpha)))
            + sqrt_6
            * (18 + alpha * (-1139 + alpha * (-379 + 75 * alpha * (35 + 18 * alpha))))
        )
        + cr
        * V2**2
        * (
            -420 * sqrt_alpha * sqrt_1_2alpha
            - 480 * V2 * sqrt_alpha * sqrt_1_2alpha
            + 9330 * alpha * sqrt_alpha * sqrt_1_2alpha
            - 1895 * V2 * alpha * sqrt_alpha * sqrt_1_2alpha
            + 16770 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            + 2385 * V2 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            + 750 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            + 16250 * V2 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            - 4500 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            + 15000 * V2 * alpha**4 * sqrt_alpha * sqrt_1_2alpha
            - sqrt_6
            * V2
            * (1 + 3 * alpha)
            * (-66 + alpha * (637 + 50 * alpha * (71 + 66 * alpha)))
            + sqrt_6
            * (-18 + alpha * (1589 + alpha * (2329 - 75 * alpha * (19 + 42 * alpha))))
        )
        + cr**3
        * (
            -240 * sqrt_alpha * sqrt_1_2alpha
            - 60 * V2 * sqrt_alpha * sqrt_1_2alpha
            + 6870 * alpha * sqrt_alpha * sqrt_1_2alpha
            + 5985 * V2 * alpha * sqrt_alpha * sqrt_1_2alpha
            + 9270 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            + 17745 * V2 * alpha**2 * sqrt_alpha * sqrt_1_2alpha
            - 2250 * V2 * alpha**3 * sqrt_alpha * sqrt_1_2alpha
            - sqrt_6
            * (
                6
                + alpha * (-163 + 18 * alpha * (49 + 100 * alpha))
                + V2
                * (1 + 3 * alpha)
                * (-54 + alpha * (299 + 25 * alpha * (103 + 54 * alpha)))
            )
        )
    )


###############################
# The coefficient \BB_5.
def b5b5(b, api):
    alpha = api.shared(b**2, label="alpha")
    cr, V2 = make_cr_and_V2(alpha, api)
    sqrt_alpha = b
    sqrt_6 = api.shared(api.sqrt(6, label="6"))
    sqrt_1_2alpha = api.shared(api.sqrt(1 + 2 * alpha, label="1 + 2 alpha"))

    return -(
        cr**4
        * sqrt_alpha
        * (1 + 3 * alpha)
        * (
            5 * sqrt_6
            + 30 * sqrt_6 * alpha
            - 168 * sqrt_alpha * sqrt_1_2alpha
        )
        + cr**2
        * V2
        * (
            15 * sqrt_6 * (-4 + 5 * V2) * sqrt_alpha
            + 5 * sqrt_6 * (-1 + 123 * V2) * alpha * sqrt_alpha
            + 10 * sqrt_6 * (14 + 153 * V2) * alpha**2 * sqrt_alpha
            + 60 * sqrt_6 * (1 + 18 * V2) * alpha**3 * sqrt_alpha
            - 36 * (-1 + V2) * sqrt_1_2alpha
            - 2 * (408 + 191 * V2) * alpha * sqrt_1_2alpha
            - (1497 + 572 * V2) * alpha**2 * sqrt_1_2alpha
            + 150 * (-3 + 7 * V2) * alpha**3 * sqrt_1_2alpha
            + 900 * V2 * alpha**4 * sqrt_1_2alpha
        )
        + cr
        * V2**2
        * (
            5 * sqrt_6 * (18 - 19 * V2) * sqrt_alpha
            + 5 * sqrt_6 * (25 - 149 * V2) * alpha * sqrt_alpha
            - 10 * sqrt_6 * (11 + 188 * V2) * alpha**2 * sqrt_alpha
            - 60 * sqrt_6 * (4 + 25 * V2) * alpha**3 * sqrt_alpha
            + 36 * (-1 + V2) * sqrt_1_2alpha
            + 7 * (71 + 45 * V2) * alpha * sqrt_1_2alpha
            + (1040 + 1021 * V2) * alpha**2 * sqrt_1_2alpha
            + 50 * (1 + 34 * V2) * alpha**3 * sqrt_1_2alpha
            + 300 * (-1 + 5 * V2) * alpha**4 * sqrt_1_2alpha
        )
        + V2**3
        * (
            40 * sqrt_6 * (-1 + V2) * sqrt_alpha
            + 75 * sqrt_6 * (-1 + 4 * V2) * alpha * sqrt_alpha
            + 20 * sqrt_6 * (3 + 37 * V2) * alpha**2 * sqrt_alpha
            + 60 * sqrt_6 * (3 + 10 * V2) * alpha**3 * sqrt_alpha
            - 12 * (-1 + V2) * sqrt_1_2alpha
            - 2 * (13 + 86 * V2) * alpha * sqrt_1_2alpha
            - (111 + 908 * V2) * alpha**2 * sqrt_1_2alpha
            - 100 * (1 + 21 * V2) * alpha**3 * sqrt_1_2alpha
            - 300 * (1 + 6 * V2) * alpha**4 * sqrt_1_2alpha
        )
        + cr**3
        * (
            -12 * sqrt_1_2alpha
            + 12 * V2 * sqrt_1_2alpha
            + 345 * alpha * sqrt_1_2alpha
            + 407 * V2 * alpha * sqrt_1_2alpha
            + 468 * alpha**2 * sqrt_1_2alpha
            + 1063 * V2 * alpha**2 * sqrt_1_2alpha
            - 150 * V2 * alpha**3 * sqrt_1_2alpha
            - 5
            * sqrt_6
            * sqrt_alpha
            * (
                -2
                + 9 * alpha
                + 18 * alpha**2
                + V2
                * (1 + 3 * alpha)
                * (5 + 4 * alpha * (7 + 3 * alpha))
            )
        )
    )

##############################################
# Build the six verification problems from the explicit functions
# and the corresponding lower/upper bounds in bounds.txt.
def _build_bernstein_problems() -> List[BoundProblem]:
    api = ProblemAPI(default_interval=B_INTERVAL)
    b = api.b
    # Each entry says: internal name, printed LaTeX symbol, formula builder,
    # lower-bound label, upper-bound label.
    entries = [
        ("b5b0", r"\tilde{\BB}_0", b5b0, "BB0_lower", "BB0_upper"),
        ("b5b1", r"\tilde{\BB}_1", b5b1, "BB1_lower", "BB1_upper"),
        ("b5b2", r"\tilde{\BB}_2", b5b2, "BB2_lower", "BB2_upper"),
        ("b5b3", r"\tilde{\BB}_3", b5b3, "BB3_lower", "BB3_upper"),
        ("b5b4", r"\tilde{\BB}_4", b5b4, "BB4_lower", "BB4_upper"),
        ("b5b5", r"\tilde{\BB}_5", b5b5, "BB5_lower", "BB5_upper"),
    ]
    return [
        BoundProblem(
            name=name,
            symbol=symbol,
            expr=formula(b, api),
            interval=tuple(B_INTERVAL),
            cuts=tuple(CUTS),
            lower_label=lower_label,
            upper_label=upper_label,
        )
        for name, symbol, formula, lower_label, upper_label in entries
    ]

##############################################
# Routine called by script.py.  It verifies the bounds for all six Bernstein coefficients
#  and returns a success flag with the list of constants used.
def bernstein_bounds():  # (Lemma \ref{lemma:Bernstein:q} and Table \ref{table:bb})
    """Verify the six two-sided bounds displayed in Table~\\ref{table:bb}."""
    lemma_label = r"Lemma~\ref{lemma:Bernstein:q} / Table~\ref{table:bb}"
    problems = _build_bernstein_problems()
    bound_labels: List[str] = []
    verified = True

    for idx, problem in enumerate(problems, start=1):
        lower_label = problem.lower_label
        upper_label = problem.upper_label
        bound_labels.extend([lower_label, upper_label])
        lower_bound = BOUNDS[lower_label]
        upper_bound = BOUNDS[upper_label]
        # Prove the current coefficient bounds.
        result = prove_table_row(
            problem,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            progress_every=(parameters.VERBOSE_COUNTER_I if parameters.VERBOSE else PROGRESS_EVERY),
            print_subintervals=(parameters.PRINT_SUBINTERVALS or parameters.VERBOSE >= 2),
        )
        row_ok = bool(result["verified"])
        verified = verified and row_ok

        if row_ok:
            lo = result["global_lower_bound"]
            hi = result["global_upper_bound"]
            assert o(lo)
            assert o(hi)
            print_iter_N(
                idx,
                total=len(problems),
                text=(
                    f"{problem.name}: {BOUNDS_STR[lower_label]} < "
                    f"{problem.symbol} < {BOUNDS_STR[upper_label]} "
                    f"(boxes={result['box_count']}, depth={result['max_depth_reached']})"
                ),
            )
        else:
            print_iter_N(
                idx,
                total=len(problems),
                text=f"{problem.name}: FAILED; first failures={result['failures'][:3]}",
            )

    bounds = [
        (BOUNDS_STR[label], label, BOUNDS_META[label]["side"])
        for label in bound_labels
    ]
    return verified, lemma_label, bounds

##############################################
# Sanity check: every label referenced by the generated problems must exist in
# bounds.txt with the expected lower/upper side.
def substituting_estimates():
    """Check that all table-bound labels used by the proof are present."""
    problems = _build_bernstein_problems()
    missing = []
    wrong_side = []
    for problem in problems:
        for label, side in ((problem.lower_label, "lower"), (problem.upper_label, "upper")):
            if label not in BOUNDS:
                missing.append(label)
                continue
            if BOUNDS_META[label]["side"] != side:
                wrong_side.append((label, BOUNDS_META[label]["side"], side))
    if missing or wrong_side:
        return False, (
            "bound-label substitution check failed: "
            f"missing={missing}, wrong_side={wrong_side}"
        )
    return True, r"Substitutions for Table~\ref{table:bb}"
