# Location-Scale CDF Derivatives at Any Order

Closed-form derivatives of \\F\\ of any order up to four, for a family
that is location-scale in its first two parameters. The general form of
[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md),
which stops at second order.

## Usage

``` r
loc_scale_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from `distrib` whose first two parameters are a
  location and a scale, and which supplies response derivatives to the
  order asked for.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, the location first and the scale second.

- order:

  The derivative order, 1 to 4.

## Value

A named list of numeric vectors, derivatives of \\F\\ itself on the
natural scale, keyed as
\[`deriv_names(distrib@params[1:2], order)`\]\[deriv_names\].

## The construction

With \\z = (q-\mu)/\sigma\\ the distribution function is \\F(q) =
F_0(z)\\, so every derivative in \\(\mu, \sigma)\\ is one Faa di Bruno
pass over that composition. The inner derivatives are \$\$F_0^{(m)}(z) =
\sigma^{m}\\\frac{\partial^{m} F}{\partial q^{m}}, \qquad
\frac{\partial^{m} F}{\partial q^{m}} = f(q)\\B\_{m-1},\$\$ with
\\B\_{m-1}\\ the complete Bell polynomial in the response derivatives of
\\\log f\\. The outer map is \\\partial^{i+j}
z/\partial\mu^{i}\partial\sigma^{j}\\, which vanishes for \\i \ge 2\\
because \\z\\ is linear in the location, so most of its partials are
exact zeros.

## What it depends on

The response derivatives of orders 3 and 4 are what make this reach past
second order: with only
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
the construction stops where
[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
stops. At orders 1 and 2 the two agree to the last bit, which is the
check that licenses the orders above.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density, \\F_0\\ the standardized
distribution function and \\B_m\\ the complete Bell polynomial.

## See also

[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
for orders 1 and 2;
[`loc_scale_deriv_cdf_k()`](https://statmodels7.github.io/distributions7/reference/loc_scale_deriv_cdf_k.md),
which registers this as a method;
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
for the partition sum;
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md).
