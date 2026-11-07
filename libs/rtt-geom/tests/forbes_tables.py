"""Generates forbes_tables.inc: values of the Forbes polynomials Q-con and Q-bfs and of their
derivatives at x = k/64, k = 0 ... 64, each exact and rounded once to double (ADR 0043,
acceptance 2). Standard library only. Usage: python forbes_tables.py > forbes_tables.inc

Independent of the recurrences that rtt-geom uses (Forbes 2010a, Eq. (5.2); Forbes 2010b,
Eqs. (A.14)-(A.16)): both families come from Gram-Schmidt over their defining inner products.

- Q-con (G. W. Forbes, Opt. Express 15(8), 5218-5226 (2007), Eqs. (5), (7), (8) and the text
  below Eq. (9)): orthogonal on [0, 1] with the weight x^4, Q_m of degree m, scaled to the
  maximum value 1 at x = 1. With <x^i, x^j> = 1/(i + j + 5) everything is an exact fraction.
- Q-bfs (G. W. Forbes, Opt. Express 18(19), 19700-19712 (2010), Eqs. (2.2)-(2.4); Forbes 2007,
  Eqs. (12), (13)): orthonormal in the slope, (2/pi) int_0^1 S_m S_n / sqrt(1 - u^2) du = d_mn
  with S_m(u) = d/du [u^2 (1 - u^2) Q_m(u^2)], and Q_m(0) > 0 (Forbes 2007, Eq. (14)). S_m S_n
  is even in u, and (2/pi) int_0^1 u^(2k) / sqrt(1 - u^2) du = (2k - 1)!! / (2k)!!, so the Gram
  matrix is exact; the normalization takes one square root, done with 50 decimal digits.

Content guards (assertions): Forbes 2007, Eq. (9) (Q-con 0 ... 5) and Eq. (14) (Q-bfs 0 ... 5)
as exact polynomials.
"""

from __future__ import annotations

import io
import sys
from collections.abc import Callable
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

getcontext().prec = 50
QCON_MAX = 20
QBFS_MAX = 16
XS = [Fr(k, 64) for k in range(65)]

Poly = list[Fr]  # coefficients in x, index = power


def evaluate(p: Poly, x: Fr) -> Fr:
    acc = Fr(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def derivative(p: Poly) -> Poly:
    return [i * p[i] for i in range(1, len(p))] or [Fr(0)]


def add(p: Poly, q: Poly, scale: Fr) -> Poly:
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else Fr(0)) + scale * (q[i] if i < len(q) else Fr(0))
            for i in range(n)]


def gram_schmidt(m_max: int, inner: Callable[[Poly, Poly], Fr]) -> list[Poly]:
    """Orthogonal (not normalized) polynomials of degree 0 ... m_max from the monomials."""
    out: list[Poly] = []
    for m in range(m_max + 1):
        p: Poly = [Fr(0)] * m + [Fr(1)]
        for q in out:
            p = add(p, q, -inner(p, q) / inner(q, q))
        out.append(p)
    return out


def qcon_inner(p: Poly, q: Poly) -> Fr:
    return sum((a * b / (i + j + 5) for i, a in enumerate(p) for j, b in enumerate(q)), Fr(0))


def slope_poly(q: Poly) -> Poly:
    """S(u) = d/du [u^2 (1 - u^2) q(u^2)] as coefficients in u."""
    s = [Fr(0)] * (2 * len(q) + 4)
    for i, c in enumerate(q):  # u^2 (1 - u^2) x^i = u^(2i+2) - u^(2i+4)
        s[2 * i + 2] += c
        s[2 * i + 4] -= c
    return [k * s[k] for k in range(1, len(s))]


def chebyshev_moment(k: int) -> Fr:
    """(2/pi) int_0^1 u^(2k) / sqrt(1 - u^2) du = (2k - 1)!! / (2k)!!."""
    out = Fr(1)
    for j in range(1, k + 1):
        out *= Fr(2 * j - 1, 2 * j)
    return out


def qbfs_inner(p: Poly, q: Poly) -> Fr:
    sp, sq = slope_poly(p), slope_poly(q)
    total = Fr(0)
    for i, a in enumerate(sp):
        for j, b in enumerate(sq):
            if a and b:
                assert (i + j) % 2 == 0  # odd times odd
                total += a * b * chebyshev_moment((i + j) // 2)
    return total


def qcon() -> list[Poly]:
    out = []
    for p in gram_schmidt(QCON_MAX, qcon_inner):
        scale = 1 / evaluate(p, Fr(1))
        out.append([c * scale for c in p])
    return out


def qbfs() -> list[tuple[Poly, Decimal]]:
    """(q, 1/sqrt(<q, q>)) with q(0) > 0: Q_m = q * the factor."""
    out = []
    for p in gram_schmidt(QBFS_MAX, qbfs_inner):
        if p[0] < 0:
            p = [-c for c in p]
        norm = qbfs_inner(p, p)
        out.append((p, 1 / (Decimal(norm.numerator) / Decimal(norm.denominator)).sqrt()))
    return out


def guards(con: list[Poly], bfs: list[tuple[Poly, Decimal]]) -> None:
    def poly(*c: int) -> Poly:
        return [Fr(v) for v in c]

    # Forbes 2007, Eq. (9).
    eq9 = [poly(1), poly(-5, 6), poly(15, -42, 28), poly(-35, 168, -252, 120),
           poly(70, -504, 1260, -1320, 495), poly(-126, 1260, -4620, 7920, -6435, 2002)]
    for m, p in enumerate(eq9):
        assert con[m] == p, f"Q-con {m}"
    # Forbes 2007, Eq. (14): q times a square root; compare q / q(0) and the squared factor.
    eq14 = [(poly(1), Fr(1)), (poly(13, -16), Fr(1, 19)),
            (poly(29, -100, 76), Fr(2, 95)),
            (poly(207, -1260, 2308, -1280), Fr(2, 2545)),
            (poly(7737, -74448, 236192, -299008, 130304), Fr(1, 9 * 131831)),
            (poly(66657, -906816, 4330400, -9186304, 8873216, -3182592), Fr(1, 9 * 6632213))]
    for m, (p, square) in enumerate(eq14):
        q, factor = bfs[m]
        ratio = p[0] / q[0]
        assert [c * ratio for c in q] == p, f"Q-bfs {m}"
        f = factor / (Decimal(ratio.numerator) / Decimal(ratio.denominator))
        s = Decimal(square.numerator) / Decimal(square.denominator)
        assert abs(f * f - s) < Decimal(10) ** -45, f"Q-bfs {m} factor"


def table(rows: list[list[float]]) -> str:
    lines = []
    for row in rows:
        lines.append("    {" + ", ".join(repr(v) for v in row) + "},")
    return "\n".join(lines)


def main() -> None:
    con = qcon()
    bfs = qbfs()
    guards(con, bfs)
    con_value = [[float(evaluate(p, x)) for x in XS] for p in con]
    con_slope = [[float(evaluate(derivative(p), x)) for x in XS] for p in con]

    def dec(v: Fr) -> Decimal:
        return Decimal(v.numerator) / Decimal(v.denominator)

    bfs_value = [[float(dec(evaluate(p, x)) * f) for x in XS] for p, f in bfs]
    bfs_slope = [[float(dec(evaluate(derivative(p), x)) * f) for x in XS] for p, f in bfs]
    out = sys.stdout
    if isinstance(out, io.TextIOWrapper):
        out.reconfigure(newline="\n")  # LF also on Windows
    out.write("// Generated by forbes_tables.py (ADR 0043, acceptance 2); do not edit.\n")
    out.write("// Q-con and Q-bfs and their derivatives d/dx at x = k/64, k = 0 ... 64, each exact\n")
    out.write("// and rounded once to double.\n\n")
    for name, rows in (("kQConValue", con_value), ("kQConDerivative", con_slope),
                       ("kQBfsValue", bfs_value), ("kQBfsDerivative", bfs_slope)):
        out.write(f"inline constexpr double {name}[{len(rows)}][65] = {{\n{table(rows)}\n}};\n\n")


if __name__ == "__main__":
    main()
