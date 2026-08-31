# Polynomial Multiplication — Karatsuba

> **Reference solution by the FHERMA team** · answers
> [`polynomial-multiplication` / `negacyclic@1.0.0`](https://fherma.io/kernels/polynomial-multiplication/specifications/negacyclic)
> on the FHERMA kernel catalogue.

`c = a·b` in the negacyclic ring `Z_q[X]/(X^N + 1)` with Karatsuba's recursion:
split each operand in half and form the product from **three** half-sized
products instead of four, recursively — `O(N^log2 3) ≈ O(N^1.585)` coefficient
multiplications rather than `O(N²)`. Below a small cutoff the plain double loop
is faster, so the leaves are schoolbook; the negacyclic fold (`X^N = −1`) and one
reduction modulo `q` per coefficient are applied once at the end. The operands
are public: this is raw compute, not encrypted. A **reference baseline** — faster
than schoolbook by a constant in the exponent, but still not competitive at
`N = 65,536`, where a real entry uses a number-theoretic transform.

```text
kernel polymul<
    N: u32,                    // ring dimension (a power of two)
    W: u32,                    // coefficient width in bits (~ log2 q)
    L: u32 = (W + 31) / 32,    // limbs per coefficient — derived
    q: tensor<L x u32> (N, W), // modulus, derived from (N, W); given at setup
>(
    %a: tensor<N x L x u32>,   // input polynomials, coefficients in [0, q)
    %b: tensor<N x L x u32>,
) -> %c: tensor<N x L x u32>   // a·b in Z_q[X]/(X^N + 1)
```

## Layout

| | |
| --- | --- |
| [`solution/`](solution/) | the measured project: `solve.py` (`init` / `run` / `free`), `fherma.toml`, generated harness |
| [`solution/README.md`](solution/README.md) | what it does and how to run |

## Running it

```sh
cd solution
python main.py <point-dir>
```

The point directory is one the specification's testing bundle produces (`make`):
it holds the point, the cases and `point/q.bin`. At measurement the platform lays
its own copy of every generated file over the clone — the only authored file is
`solve.py`.

## License

MIT.
