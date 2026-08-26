# Inverse Gaussian Log-CDF Gradient in Mean and Rate

Closed form, by the chain rule on the dispersion parametrization's
gradient through \\\phi = 1/\lambda\\. The mean component is unchanged;
the rate component is the dispersion one times \\-1/\lambda^2\\.

## Arguments

- distrib:

  An `InvGauss2Distrib` object, from
  [`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

- q:

  A numeric vector of quantiles, positive.

- theta:

  A named list with components `mu` (positive) and `lambda` (positive),
  each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `lambda`, each the length
of `q` recycled against `theta`.

## Details

The inverse Gaussian is unusual among positive families in having an
elementary distribution function, \\F(y) = \Phi(a) +
e^{2/(\phi\mu)}\\\Phi(b)\\, so the parent's derivatives are closed and
the chain carries them.

The second and higher orders are registered elsewhere, by
[`register_mapped_cdf_k()`](https://statmodels7.github.io/distributions7/reference/register_mapped_cdf_k.md)
in `cdf_mapped_higher.R`, at orders 2 to 4. That file's registration
supersedes anything this one would give: the Hessian was left on the
fallback while the parent differenced its own second order, and moved
onto the chain when the parent stopped.

## Notation

\\\mu \> 0\\ is the mean, \\\lambda \> 0\\ the rate, \\\phi =
1/\lambda\\ the dispersion and \\\Phi\\ the standard normal distribution
function.

## See also

[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for the generic;
[`mapped_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv.md)
and
[`chain_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/chain_cdf_deriv.md);
[`invgauss2_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss2_distrib.md).

## Examples

``` r
d <- invgauss2_distrib()
q <- c(0.5, 2, 5)
th <- list(mu = 2, lambda = 8)

# Against a central difference of the cdf, which shares no arithmetic.
fd <- numerical_cdf_deriv(d, q, th, order = 1)
max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) / unlist(fd) - 1))
#> [1] 2.463587e-10
```
