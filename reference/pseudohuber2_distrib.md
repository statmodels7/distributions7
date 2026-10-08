# Pseudo-Huber Distribution, Standard-Deviation Parametrization

Creates a pseudo-Huber distribution object parametrized by the location
\\\mu\\, the standard deviation \\\sigma\\ and the shape \\\nu\\.

## Usage

``` r
pseudohuber2_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_nu = log_link()
)
```

## Arguments

- link_mu:

  A link for the location, by default
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html).

- link_sigma:

  A link for the standard deviation, by default
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_nu:

  A link for the shape, by default
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class `PseudoHuber2Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"pseudo huber2"`, `params`
`c("mu", "sigma", "nu")` and `link_params` the three links given here.

## Details

The family is that of
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md),
whose scale \\\sigma_1\\ is replaced by the standard deviation \\\sigma
= \sigma_1 \sqrt{R(\nu)}\\, \\R(\nu) = \sqrt\nu\\
K_2(\sqrt\nu)/K_1(\sqrt\nu)\\. As \\\nu \to \infty\\ the distribution
tends to the gaussian with standard deviation \\\sigma\\, so the
gaussian limit is reached along the single coordinate \\\nu\\ at a fixed
\\\sigma\\. In the scale parametrization the same limit needs \\\sigma_1
\to 0\\ together with \\\nu \to \infty\\, and the correlation of
\\(\log\sigma_1, \log\nu)\\ in the inverse expected information is -0.93
at \\\nu = 1\\ and -0.9997 at \\\nu = 1000\\; in this parametrization
the correlation of \\(\log\sigma, \log\nu)\\ is -0.37 and -0.05.

Every derivative of the log-density in the parameters to order five is
in closed form (see
[`distrib_gradient.PseudoHuber2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.PseudoHuber2Distrib.md)),
the expected information is computed by one quadrature per distinct
shape, and the distribution function by quadrature of the density.

## References

Barndorff-Nielsen, O. (1978). Hyperbolic distributions and distributions
on hyperbolae. *Scandinavian Journal of Statistics*, **5**(3), 151-157.

## See also

[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
for the scale parametrization;
[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md)
for the Student t parametrized by its standard deviation;
[PseudoHuber2Distrib](https://statmodels7.github.io/distributions7/reference/PseudoHuber2Distrib.md)
for the class.

## Examples

``` r
d <- pseudohuber2_distrib()
th <- list(mu = 0.4, sigma = 1.5, nu = 2)
# sigma is the standard deviation
c(variance = variance(d, th), sigma2 = 1.5^2)
#> variance   sigma2 
#>     2.25     2.25 
# a large shape is the gaussian with standard deviation sigma
yy <- c(-1, 0.5, 2)
all.equal(distrib_pdf(d, yy, list(mu = 0.4, sigma = 1.5, nu = 1e12)),
          dnorm(yy, 0.4, 1.5), tolerance = 1e-6)
#> [1] TRUE
```
