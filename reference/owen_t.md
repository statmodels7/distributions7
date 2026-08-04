# Owen's T Function

Computes \\T(h, a) = \dfrac{1}{2\pi}\displaystyle\int_0^{a}
\dfrac{e^{-h^2(1 + x^2)/2}}{1 + x^2}\\dx\\, which is what the skew
normal distribution function is written in.

## Usage

``` r
owen_t(h, a)
```

## Arguments

- h:

  A numeric vector.

- a:

  A numeric vector, recycled against `h`.

## Value

A numeric vector.

## Details

The integral is one-dimensional over a finite range with a bounded,
smooth integrand, so adaptive quadrature evaluates it to near machine
precision. That is a better object to integrate than the density over a
semi-infinite range, which is what the base class would otherwise do.

Two identities keep the extremes exact rather than quadrature-bound:
\\T(h, a) = -T(h, -a)\\, and \\T(h, \infty) =
\tfrac{1}{2}\Phi(-\|h\|)\\.

## References

Owen, D. B. (1956). Tables for computing bivariate normal probabilities.
*Annals of Mathematical Statistics* 27, 1075-1090.
