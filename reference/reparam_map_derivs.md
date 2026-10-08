# Explicit Map Derivatives of the Second Parametrizations

Each function returns, for its family's map \\\theta = h(\psi)\\, the
non-zero partials \\\partial^B \theta_i / \partial \psi_B\\ up to order
`order`: a list over parent parameters, each a list keyed by the sorted
tuple of \\\psi\\ positions. A missing key is an exact zero. A partial
of higher order than `order` is not formed. Every formula is derived by
hand and validated against one numerical pass per order in the tests.

## Usage

``` r
md_betabinom1(psi, order)

md_lognormal2(psi, order)

md_weibull3(psi, order)

md_student_t2(psi, order)

md_gengamma2(psi, order)

md_invgauss2(psi, order)

md_skewnormal2(psi, order)

md_laplace2(psi, order)

md_gaussian2(psi, order)

md_gaussian3(psi, order)
```

## Arguments

- psi:

  The aligned list of the new parameters.

- order:

  The highest order of partial to form, an integer from 1 to 4.

## Value

A list over parent parameters of keyed partial tables.
