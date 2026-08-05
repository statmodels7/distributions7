# Higher Derivatives of the Beta-Binomial in Its Shapes

The log-mass function is a sum of log-gamma terms, so its derivative of
order \\k\\ is the same sum with \\\psi^{(k-1)}\\ in place of
\\\log\Gamma\\. All four orders follow from one routine.

## Usage

``` r
betabinom2_component(y, a, b, n, k, i, j)
```

## Arguments

- y:

  A numeric vector of observations.

- a, b:

  The two shapes.

- n:

  The size.

- k:

  The polygamma order, `order - 1`.

- i:

  The number of \\a\\ indices in the component.

- j:

  The number of \\b\\ indices.

## Value

A numeric vector.

## Details

With \\a\\ and \\b\\ the shapes and \\n\\ the size, the log-mass is
\\\log\Gamma(y+a) + \log\Gamma(n-y+b) - \log\Gamma(n+a+b) -
\log\Gamma(a) - \log\Gamma(b) + \log\Gamma(a+b)\\ up to a constant. A
derivative in \\a\\ alone differentiates the first, third, fourth and
sixth terms; one in \\b\\ alone the second, third, fifth and sixth; a
mixed one only the two terms carrying \\a+b\\.

## See also

[`betabinom2_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
