# Third and Fourth Derivatives of a Compiled Distribution Function

The derivatives of \\F\\ of order three or four in the parameters, for
the families whose compiled distribution function integrates its shape
derivatives
([`compiled_cdf()`](https://statmodels7.github.io/distributions7/reference/compiled_cdf.md)):
gamma1, gamma2, chisq, both generalized gammas, both betas, both von
Mises, both Student t, both pseudo-Huber and the skew t.

## Usage

``` r
compiled_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  A family named in the description.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters on the parameter scale.

- order:

  The derivative order, 3 or 4.

## Value

A named list of numeric vectors, derivatives of \\F\\ on the natural
scale and the lower tail, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md).

## Details

With \\f = e^{\ell}\\, the derivative of \\f\\ in a multi-index \\I\\ is
\\f B_I\\, where \\B_I\\ is the complete Bell polynomial in the
derivatives of \\\ell\\
([`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)),
so \\\partial_I F(q) = \int\_{lo}^{q} f B_I\\, taken as
\\-\int\_{q}^{hi} f B_I\\ on the side `cdf_rule_cpp()` chooses, by the
nodes and weights of the rule of the first two orders. The family's own
third and fourth derivatives of \\\ell\\ supply \\B_I\\.

A family location-scale in its first two parameters (both Student t,
both pseudo-Huber, the skew t) reads its components in the location and
the scale at \\q\\: with \\\partial\_\mu F = -f(q)\\ and
\\\partial\_\sigma F = -(q - \mu) f(q)/\sigma\\, \\\partial\_{\mu J} F =
-f(q) B_J\\ and, for a multi-index with \\m\\ scale indices and shape
indices \\S\\, \\\sum\_{j=0}^{m-1} \binom{m-1}{j} c_j f(q)
B\_{\sigma^{m-1-j} S}\\ with \\c_j = -(q - \mu)(-1)^j j!/\sigma^{j+1}\\.

## See also

[`numerical_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/numerical_cdf_deriv_k.md),
the stencil this replaces for these families;
[`continuous_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/continuous_cdf_deriv_k.md),
which chooses between them.
