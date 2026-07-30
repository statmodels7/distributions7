# CDF Derivatives of a Lattice Distribution, Exactly

Evaluates \\d^I F(q)\\ for a discrete distribution as the finite sum
\\\sum\_{y \le q} f(y) \\ (d^I f / f)(y)\\.

## Usage

``` r
discrete_cdf_deriv(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"discrete_distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 or 2.

## Value

A named list of derivative component vectors of \\F\\.

## Details

This is the governing identity \\d^I F(q) / F(q) = \mathbb{E}\[d^I f / f
\mid Y \le q\]\\ written out. For a lattice family the conditional
expectation is a finite sum whenever the support has a finite lower
bound – which the discrete class requires – so the identity is exact
rather than an approximation, and it is used directly.

Note for anyone writing a test: checking this against the
partial-expectation sum proves nothing, because it is the same sum
computed twice. A discrete implementation has to be checked against
finite differences of the cdf.

## See also

[`cdf_tail_scale`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md)
