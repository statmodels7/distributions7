# Inverse Gaussian Distribution, Mean and Shape

Builds the distribution object for the inverse Gaussian family on \\(0,
\infty)\\ in its classical parametrization, the mean \\\mu \> 0\\ and
the shape \\\lambda \> 0\\, so that \\\operatorname{Var}(Y) =
\mu^3/\lambda\\. The returned object carries closed-form derivatives of
the log-density to fourth order, in the parameters and in the response,
and closed-form moments, so every generic of the toolkit answers without
a numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both parameters being positive.

## Usage

``` r
invgauss2_distrib(link_mu = log_link(), link_lambda = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  mean positive.

- link_lambda:

  A `link` object from `linkfunctions7` for the shape \\\lambda\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `InvGauss2Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"invgauss2"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "lambda")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## The parametrization

The density on \\y \in (0, \infty)\\ is \$\$f(y; \mu, \lambda) =
\sqrt{\dfrac{\lambda}{2\pi y^{3}}}
\exp\left\\-\dfrac{\lambda(y-\mu)^{2}}{2\mu^{2}y}\right\\,\$\$ the mean
is \\\mu\\, the variance \\\mu^3/\lambda\\, the skewness
\\3\sqrt{\mu/\lambda}\\ and the excess kurtosis \\15\mu/\lambda\\. The
distribution function is elementary in the standard normal one and the
quantile function is its numerical inverse.

This is the same law as
[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md),
which carries a dispersion \\\phi = 1/\lambda\\. The map between the two
moves one coordinate at a time, the mean being untouched, so both sets
of derivatives stay elementary. The shape is the parametrization the
family is usually written in, and it is the one in which the law is the
first-passage time of a Brownian motion with drift \\\mu\\ and variance
parameter \\\lambda\\.

## Derivatives

The score is \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{\lambda(y-\mu)}{\mu^3}, \qquad \dfrac{\partial \ell}{\partial
\lambda} = \dfrac{1}{2\lambda} - \dfrac{(y-\mu)^2}{2\mu^2 y},\$\$ and
the expected Hessian is \$\$\mathbb{E}\left\[\ell^{(\mu\mu)}\right\] =
-\dfrac{\lambda}{\mu^3}, \qquad
\mathbb{E}\left\[\ell^{(\mu\lambda)}\right\] = 0, \qquad
\mathbb{E}\left\[\ell^{(\lambda\lambda)}\right\] =
-\dfrac{1}{2\lambda^2}.\$\$ The zero off-diagonal makes the mean and the
shape orthogonal.

The log-density is **linear in \\\lambda\\** apart from
\\\tfrac12\log\lambda\\, and that shapes every higher order. The pure
shape derivatives are those of the logarithm alone, \\-1/(2\lambda^2)\\,
\\1/\lambda^3\\ and \\-3/\lambda^4\\; every component naming the shape
twice or more alongside the mean is exactly zero; and each remaining
component is at most linear in the response, so every expectation needs
only \\\mathbb{E}\[Y\] = \mu\\ and all four orders are closed form. It
also makes the pure shape entry of the observed Hessian a negative
constant, where the corresponding entry in the dispersion
parametrization can be positive.

Third and fourth orders are in
[`distrib_deriv3.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGauss2Distrib.md)
and
[`distrib_deriv4.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.InvGauss2Distrib.md).
The derivatives in the response and the mixed block are registered
elsewhere in the package and are closed form as well; see
[`distrib_grad_y.InvGauss2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.InvGauss2Distrib.md).

## Estimation

Both maximum likelihood estimates are available in closed form:
\$\$\hat\mu = \bar y, \qquad \dfrac{1}{\hat\lambda} = \dfrac{1}{n}\sum_i
\dfrac{1}{y_i} - \dfrac{1}{\bar y}.\$\$
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches them numerically on the link scale, and the example below checks
both against the sample.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu \> 0\\ the mean
and \\\lambda \> 0\\ the shape, with \\\operatorname{Var}(Y) =
\mu^3/\lambda\\. Here \\\lambda\\ is this family's shape and not a
penalty parameter or an eigenvalue. \\\eta\\ is a parameter on the
unconstrained scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Chhikara, R. S. and Folks, J. L. (1989). *The Inverse Gaussian
Distribution: Theory, Methodology, and Applications*. Marcel Dekker, New
York.

## See also

[`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
for the same law in the mean and a dispersion, which is the generalized
linear model parametrization;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
and
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md)
for other positive families;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for the counts this law mixes a Poisson into;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[InvGauss2Distrib](https://statmodels7.github.io/distributions7/reference/InvGauss2Distrib.md)
for the class.

## Examples

``` r
d <- invgauss2_distrib()
d
#> Distribution: Invgauss2
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu     (mean)               | Link: log        | Domain: (0, Inf)
#>   lambda (shape)              | Link: log        | Domain: (0, Inf)

# The same law as invgauss1 at phi = 1/lambda.
y <- c(1, 2, 3)
th <- list(mu = 2, lambda = 3)
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(invgauss1_distrib(), y, list(mu = 2, phi = 1 / 3)))
#> [1] TRUE

# Moments: variance mu^3/lambda, skewness 3 sqrt(mu/lambda).
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>      mean       var      skew      kurt 
#>  2.000000  2.666667  2.449490 10.000000 
c(2^3 / 3, 3 * sqrt(2 / 3), 15 * 2 / 3)
#> [1]  2.666667  2.449490 10.000000

# Fitting recovers the closed-form maximum likelihood estimates.
set.seed(3)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted = coef(fit),
      closed = c(mu = mean(z),
                 lambda = 1 / (mean(1 / z) - 1 / mean(z))))
#>              mu   lambda
#> fitted 1.959755 3.018427
#> closed 1.959755 3.018428

# The mean and the shape are orthogonal: the mixed entry is 0.
distrib_expected_hessian(d, 1, th)$mu_lambda
#> [1] 0
```
