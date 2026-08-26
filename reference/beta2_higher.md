# Higher Derivatives of the Beta in Its Shapes

Returns the derivatives of the beta log-density in the two shapes at
order 2, 3 or 4. Each component is a difference of polygamma functions
at \\\alpha\\, \\\beta\\ and \\\alpha+\beta\\: with \\k = \\`order - 1`,
a component naming only \\\alpha\\ is \\\psi^{(k)}(\alpha+\beta) -
\psi^{(k)}(\alpha)\\, one naming only \\\beta\\ is
\\\psi^{(k)}(\alpha+\beta) - \psi^{(k)}(\beta)\\, and any mixed
component is \\\psi^{(k)}(\alpha+\beta)\\.

## Usage

``` r
beta2_higher(theta, n, order)
```

## Arguments

- theta:

  A named list with components `alpha` and `beta`, each a numeric
  vector. Both must be strictly positive.

- n:

  The length to recycle each component to, normally `length(y)`.

- order:

  The derivative order, `2L`, `3L` or `4L`. Any other value falls
  through to the fourth-order branch.

## Value

A named list of numeric vectors, each of length `n`: three components at
order 2 (`alpha_alpha`, `alpha_beta`, `beta_beta`), four at order 3 and
five at order 4, named for the distinct multi-indices.

## Details

The data enter the log-density only through \\(\alpha-1)\log y +
(\beta-1)\log(1-y)\\, which is linear in the two parameters, so the
second derivative already kills it. Every order from two upwards is
therefore free of the response, the observed and the expected
derivatives are the same numbers, and no expectation is computed
anywhere in this family. The three polygamma values are evaluated once
each and recycled.

## See also

[`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
for the family, and
[`distrib_hessian.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta2Distrib.md),
[`distrib_deriv3.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Beta2Distrib.md)
and
[`distrib_deriv4.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Beta2Distrib.md)
for the methods that call it.
