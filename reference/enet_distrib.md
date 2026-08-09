# Construct an Elastic-Net Distribution

The density whose negative logarithm is the elastic-net penalty: the
product of a Laplace and a Gaussian at the same location, normalized. It
exists so that `elasticnet_penalty` is the same construction as ridge
and lasso – a `distrib_penalty` over a
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.md)
family – rather than a branch of its own.

## Usage

``` r
enet_distrib(
  link_mu = identity_link(),
  link_lambda = log_link(),
  link_alpha = logit_link()
)
```

## Arguments

- link_mu:

  A linkfunctions7 link for the location.

- link_lambda:

  A link for the overall rate, positive.

- link_alpha:

  A link for the mixing weight, in \\(0,1)\\.

## Value

An object of class
[`EnetDistrib`](https://statmodels7.github.io/distributions7/reference/EnetDistrib.md).

## Details

Writing \\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\, the
normalizing constant is \\Z = 2M(a/\sqrt{c})/\sqrt{c}\\ with \\M\\ the
Mills ratio, which is finite at both ends: \\2/a\\ as \\\alpha \to 1\\
and \\\sqrt{2\pi/c}\\ as \\\alpha \to 0\\. Every derivative in the two
rates is a polynomial in \\x = a/\sqrt{c}\\ and \\G = \mathrm{d}\log
M/\mathrm{d}x\\, and the chain to \\(\lambda, \alpha)\\ is bilinear.

Orders three and four in the parameters come from the numerical
fallback; the response derivatives are exact at every order, the third
and beyond being zero.

## See also

[`laplace2_distrib`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md),
[`gaussian1_distrib`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)

## Examples

``` r
d <- enet_distrib()
theta <- list(mu = 0, lambda = 2, alpha = 0.5)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.1701518 0.7625676 0.1701518

# the density integrates to one, and the two ends are the two families
stats::integrate(function(u) distrib_pdf(d, u, theta), -Inf, Inf)$value
#> [1] 1
distrib_pdf(d, 0.7, list(mu = 0, lambda = 2, alpha = 1 - 1e-10))
#> [1] 0.246597
distrib_pdf(laplace2_distrib(), 0.7, list(mu = 0, lambda = 2))
#> [1] 0.246597
```
