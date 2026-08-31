"""Karatsuba multiplication in Z_q[X] / (X^N + 1).

Three half-sized products instead of four, recursively: N^log2(3) ~ N^1.585
coefficient multiplications rather than N^2. The saving is in the count of
multiplications and it is paid for in additions, which is why it only wins once
the coefficients are wide enough for a multiplication to cost more than the
several additions that replaced it.
"""
import sys

from fherma import Inputs, Outputs, Point, pack, unpack, unpack_one

#: Below this, the double loop is quicker than the recursion that avoids it.
CUTOFF = 32


def init(p: Point):
    """Room to recurse, and the modulus, which is what the point decides.

    Depth is log2(N) frames of Python; the default limit of 1000 is comfortable
    at any interesting size, but a solution that died of its own recursion at
    the size it claims to handle would be reporting a limit of the interpreter
    as a property of the method. q is fixed for the point and read here, once.
    """
    if p.N > 1:
        sys.setrecursionlimit(max(1000, 20 * p.N.bit_length()))
    return p.N, p.L, unpack_one(p.q)


def run(state, inp: Inputs) -> Outputs:
    n, limbs, q = state

    a = unpack(inp.a)
    b = unpack(inp.b)

    product = _karatsuba(a, b)

    # X^N = -1, so every term at or above N folds back negatively. Done once at
    # the end rather than inside the recursion: the recursion multiplies
    # polynomials, and the ring is not its business.
    for degree in range(n, len(product)):
        product[degree - n] -= product[degree]

    return Outputs(c=pack([one % q for one in product[:n]], limbs, "u32"))


def _karatsuba(a: list, b: list) -> list:
    """The plain product of two polynomials of equal length, unreduced."""
    n = len(a)
    if n <= CUTOFF:
        return _schoolbook(a, b)

    half = (n + 1) // 2
    pad = [0] * (2 * half - n)

    a_low, a_high = a[:half], a[half:] + pad
    b_low, b_high = b[:half], b[half:] + pad

    low = _karatsuba(a_low, b_low)
    high = _karatsuba(a_high, b_high)

    a_sum = [x + y for x, y in zip(a_low, a_high)]
    b_sum = [x + y for x, y in zip(b_low, b_high)]
    middle = _karatsuba(a_sum, b_sum)
    for i, value in enumerate(low):
        middle[i] -= value
    for i, value in enumerate(high):
        middle[i] -= value

    product = [0] * (4 * half - 1)
    for i, value in enumerate(low):
        product[i] += value
    for i, value in enumerate(middle):
        product[i + half] += value
    for i, value in enumerate(high):
        product[i + 2 * half] += value

    return product[: 2 * n - 1]


def _schoolbook(a: list, b: list) -> list:
    product = [0] * (2 * len(a) - 1)
    for i, ai in enumerate(a):
        if not ai:
            continue
        for j, bj in enumerate(b):
            product[i + j] += ai * bj
    return product


def free(state) -> None:
    return None
