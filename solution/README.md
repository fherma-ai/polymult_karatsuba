# Polynomial Multiplication — Karatsuba

A reference implementation of [`polynomial-multiplication/negacyclic@1.0.0`](https://fherma.io/kernels/polynomial-multiplication) for the FHERMA catalogue.

## What it does

Computes `c = a · b` in the negacyclic ring `Z_q[X] / (X^N + 1)` with Karatsuba's recursion:
split each operand in half and form the product from **three** half-sized products instead of
four,

```
a·b = low + (mid − low − high)·X^half + high·X^{2·half}
```

where `low = a_lo·b_lo`, `high = a_hi·b_hi`, and `mid = (a_lo + a_hi)(b_lo + b_hi)`. That is
`O(N^log2 3) ≈ O(N^1.585)` coefficient multiplications rather than `O(N²)`. Below a small cutoff
the plain double loop is faster, so the leaves are schoolbook. After the recursion the negacyclic
fold (`X^N = -1`) and one reduction modulo `q` per coefficient are applied once.

`q` is fixed for the benchmark point and read once in `init`, before timing starts; the timed
`run` does only the multiplication.

This is a **reference baseline** — faster than schoolbook by a constant-in-the-exponent, but still
not competitive at `N = 65,536`, where a real entry uses a number-theoretic transform.

## Interface

Coefficients are `W`-bit integers in `[0, q)`, each stored as `L = ceil(W / 32)` little-endian
`u32` limbs; `a`, `b`, `c` are `tensor<N × L × u32>`. Parameters and layout are defined by the
specification.

## Run it

```
python main.py <point-directory>
```

The directory is one the specification's testing bundle produced (`make`). The answer is written
to `out/`.
