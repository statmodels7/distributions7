# Location-Scale CDF Derivatives at Any Order

Closed-form derivatives of \\F\\ of any order up to four, for a family
that is location-scale in its first two parameters.

## Usage

``` r
loc_scale_cdf_deriv_k(distrib, q, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, location first and scale second.

- order:

  The derivative order, 1 to 4.

## Value

A named list of derivative components of \\F\\.

## Details

With \\z = (q-\mu)/\sigma\\ the distribution function is \\F(q) =
F_0(z)\\, so every derivative in \\(\mu, \sigma)\\ is one Faa di Bruno
pass over that composition. The inner derivatives are \\F_0^{(m)}(z) =
\sigma^{m} \partial^{m} F/\partial q^{m}\\, and \\\partial^{m}
F/\partial q^{m} = f(q) B\_{m-1}\\, the complete Bell polynomial in the
response derivatives of \\\log f\\. The map is \\\partial^{i+j}
z/\partial\mu^{i}\partial\sigma^{j}\\, which vanishes for \\i \ge 2\\
because \\z\\ is linear in the location.

This is what the response derivatives of order three and four are for:
with only
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md)
the construction stops at second order, which is where
[`loc_scale_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
stops. At orders one and two the two agree exactly.

## See also

[`loc_scale_cdf_deriv`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md),
[`chain_assemble`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
