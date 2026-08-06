# Explicit Map Derivatives of the Second Parametrizations

Each function returns, for its family's map \\\theta = h(\psi)\\, the
non-zero partials \\\partial^B \theta_i / \partial \psi_B\\ to fourth
order: a list over parent parameters, each a list keyed by the sorted
tuple of \\\psi\\ positions. A missing key is an exact zero. Every
formula is derived by hand and validated against one numerical pass per
order in the tests.

## Usage

``` r
md_betabinom1(psi)

md_lognormal2(psi)

md_weibull3(psi)

md_student_t2(psi)

md_gengamma2(psi)

md_invgauss2(psi)

md_skewnormal2(psi)

md_gaussian2(psi)

md_gaussian3(psi)
```

## Arguments

- psi:

  The aligned list of the new parameters.

## Value

A list over parent parameters of keyed partial tables.
