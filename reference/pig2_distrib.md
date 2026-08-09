# Poisson-Inverse Gaussian Distribution in Its Orthogonal Parametrization

Creates a Poisson-inverse Gaussian distribution in the parametrization
whose parameters are orthogonal – the expected information is diagonal –
which is gamlss's `PIG2`. The mean stays \\\mu\\; the second parameter
\\\alpha\\ is exactly the argument the Bessel function of the mass
function is evaluated at, related to the dispersion of
[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
by \\\alpha = \sqrt{1 + 2\sigma\mu}/\sigma\\ (gamlss states the inverse
of the same map, \\\alpha = 1/(\sqrt{\mu^2 + \sigma_2^2} - \mu)\\ with
its own \\\sigma_2\\, which coincides with this \\\alpha\\).

## Usage

``` r
pig2_distrib(link_mu = log_link(), link_alpha = log_link())
```

## Arguments

- link_mu:

  The link for \\\mu\\; defaults to `log_link()`.

- link_alpha:

  The link for \\\alpha\\; defaults to `log_link()`.

## Value

A `Pig2Distrib` object.

## Details

Orthogonality makes the maximum likelihood estimates of \\\mu\\ and
\\\alpha\\ asymptotically independent, so Fisher scoring steps in one
parameter do not disturb the other; the cost is that \\\alpha\\ has no
moment reading of its own. Derivatives to fourth order are exact,
computed by the same compiled kernel as
[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
with \\\alpha\\ a seed variable of the jet, so the Bessel argument needs
no chain rule at all.

## The distribution

\$\$P(Y=y) =
\sqrt{\frac{2\alpha}{\pi}}\\\frac{\mu^{y}e^{1/\sigma}}{(\alpha\sigma)^{y}\\y!}\\K\_{y-1/2}(\alpha),
\qquad \sigma = \frac{1}{\sqrt{\mu^{2}+\alpha^{2}} - \mu}\$\$ on \\y \in
\\0, 1, \dots\\\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \mu +
\sigma\mu^{2}\$\$

## References

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Applied Statistics* 54(3),
507–554.

Heller, G. Z., Couturier, D.-L., and Heritier, S. R. (2019). Beyond mean
modelling: bias due to misspecification of dispersion in Poisson-inverse
Gaussian regression. *Biometrical Journal* 61(2), 333–342.

## See also

[`pig1_distrib`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)

## Examples

``` r
d <- pig2_distrib()
theta <- list(mu = 3, alpha = 1.2)
distrib_pdf(d, 0:5, theta)
#> [1] 0.37949984 0.21925497 0.11611779 0.06810319 0.04421579 0.03093976
mean(d, theta)
#> [1] 3
```
