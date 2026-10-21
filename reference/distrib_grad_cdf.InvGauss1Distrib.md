# Inverse Gaussian Log-CDF Derivatives

Returns the derivatives in closed form at every order from one to four.
The inverse Gaussian's distribution function is elementary, \$\$F(q) =
\Phi(a) + e^{c}\\\Phi(b), \qquad a = \frac{q/\mu - 1}{\sqrt{\phi q}},
\quad b = -\frac{q/\mu + 1}{\sqrt{\phi q}}, \quad c =
\frac{2}{\phi\mu},\$\$ so it can simply be differentiated. This page's
four registrations are made together by
[`register_phi_terms_cdf()`](https://statmodels7.github.io/distributions7/reference/register_phi_terms_cdf.md).

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- q:

  A numeric vector of quantiles, positive.

- theta:

  A named list with components `mu` (positive) and `phi` (positive);
  each is a numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors of the requested order, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md):
two components for the gradient, three for the Hessian, four at order 3
and five at order 4.

## Separable partial derivatives

Each of \\a\\, \\b\\ and \\c\\ is a product of a function of the mean
and a function of the dispersion, so every mixed partial is a product of
two one-variable derivatives. The two terms are then a Leibniz split
between the weight and the tail, which
[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
runs on the package's own partition enumeration.

## Overflow of the weight

\\e^{c}\\ overflows to `Inf` once \\c = 2/(\phi\mu)\\ exceeds about 709;
the exponent is 2000 at \\\mu = 0.01\\, \\\phi = 0.1\\. Where it
overflows, \\\Phi(b)\\ underflows, and the product is finite. The weight
and the tail are therefore combined as
`exp(c + pnorm(b, log.p = TRUE))`.

## Accuracy

The derivatives are closed forms. A product stencil on the same
distribution function loses accuracy as the order grows, so the family
registers all four orders.

## Notation

\\\mu \> 0\\ is the mean, \\\phi \> 0\\ the dispersion, \\\Phi\\ the
standard normal distribution function and \\F\\ the inverse Gaussian's.

## See also

[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
for the construction;
[`register_phi_terms_cdf()`](https://statmodels7.github.io/distributions7/reference/register_phi_terms_cdf.md),
which makes the four registrations;
[`distrib_grad_cdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.EnetDistrib.md),
the other family of this shape;
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

## Examples

``` r
d <- invgauss1_distrib()
q <- c(0.5, 2, 5)
th <- list(mu = 2, phi = 0.5)

# The gradient is compared with a central difference of the cdf, which
# shares no arithmetic.
fd <- numerical_cdf_deriv(d, q, th, order = 1)
max(abs(unlist(distrib_grad_cdf(d, q, th, log = FALSE)) / unlist(fd) - 1))
#> [1] 9.956747e-11

# The gradient is finite although the weight overflows: exp(2 / (phi mu))
# is Inf here.
exp(2 / (0.01 * 0.1))
#> [1] Inf
distrib_grad_cdf(d, 0.02, list(mu = 0.01, phi = 0.1), log = FALSE)
#> $mu
#> [1] -3.174073e-106
#> 
#> $phi
#> [1] -7.929894e-108
#> 
```
