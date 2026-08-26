# Negative Binomial Distribution, NB1

Builds the distribution object for the negative binomial family on the
non-negative integers whose variance is **linear** in the mean: with a
mean \\\mu \> 0\\ and a dispersion \\\theta \> 0\\,
\\\operatorname{Var}(Y) = \mu(1+\theta)\\. The returned object carries
closed-form derivatives of the log-mass to fourth order and closed-form
moments.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both parameters being positive.

## Usage

``` r
negbin1_distrib(link_mu = log_link(), link_theta = log_link())
```

## Arguments

- link_mu:

  A `link` object from `linkfunctions7` for the mean \\\mu\\. Defaults
  to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  mean positive.

- link_theta:

  A `link` object from `linkfunctions7` for the dispersion \\\theta\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `NegBin1Distrib`, inheriting from
`discrete_distrib`, with `distrib_name` `"negbin1"`, `dimension`
`"univariate"`, `bounds` `c(0, Inf)`, `params` `c("mu", "theta")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## Two negative binomials, and they are two families

Two negative binomials are in common use and they are **different
families**, not two parametrizations of one. Here the variance is
\\\mu(1+\theta)\\, growing in proportion to the mean, so the dispersion
relative to a Poisson is the same at every mean;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
has \\\mu + \mu^2/\theta\\, growing quadratically. Fitting one is not
fitting the other, and a likelihood ratio between them is not a test of
nested models.

The difference is visible in where the mean sits. The size is \\r =
\mu/\theta\\ and the success probability \\1/(1+\theta)\\, so \\\mu\\
appears **inside** the gamma functions of the mass; in the quadratic
form the size is \\\theta\\ and the mean stays outside them.

## The parametrization

The mass on \\y = 0, 1, 2, \ldots\\ is \$\$P(Y = y) = \dfrac{\Gamma(y +
\mu/\theta)}{\Gamma(\mu/\theta)\\y!}
\left(\dfrac{1}{1+\theta}\right)^{\mu/\theta}
\left(\dfrac{\theta}{1+\theta}\right)^{y},\$\$ the distribution function
is the partial sum and the quantile function its generalized inverse.
The mean is \\\mu\\ and the variance \\\mu(1+\theta)\\. Two settings are
worth recognizing: at \\\mu = \theta\\ the size is 1 and the law is the
geometric, and as \\\theta\\ goes to zero it tends to the Poisson, where
the quadratic form needs \\\theta \to \infty\\ instead.

## Derivatives

Everything follows from the chain rule through the size. With \\P =
\psi(y+r) - \psi(r) - \log(1+\theta)\\ and \\\psi\\ the digamma
function, the score is \$\$\dfrac{\partial \ell}{\partial \mu} =
\dfrac{P}{\theta}, \qquad \dfrac{\partial \ell}{\partial \theta} =
-\dfrac{\mu}{\theta^2}P - \dfrac{r}{1+\theta} + \dfrac{y}{\theta} -
\dfrac{y}{1+\theta},\$\$ and the Hessian is the same chain rule at
second order. In the expected information every term carrying \\P\\
drops out, its expectation vanishing by the first Bartlett identity, and
only \\\mathbb{E}\[\psi'(Y+r)\]\\ remains, which is summed against the
exact mass out to a far-tail quantile.

**The mean and the dispersion are not orthogonal here.** The mixed entry
of the expected information is small but non-zero at every setting
measured, where
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
has exactly zero there.

## What happens as theta goes to zero

The family tends to the Poisson, so every derivative in \\\theta\\ tends
to a finite limit while the pieces it is assembled from run away: the
size \\r = \mu/\theta\\ grows without bound and the chain rule divides
by powers of \\\theta\\. The score's digamma difference is computed in a
form that performs its own cancellation symbolically, and the value
converges onto \\\\(y-\mu)^2 - y\\/(2\mu)\\, holding to about five
significant figures at \\\theta = 10^{-8}\\. Two quantities are **not**
rewritten, and both pages say so:

- the expected information in \\\theta\\, which turns positive from
  about \\\theta = 10^{-6}\\ and leaves the matrix indefinite
  ([`distrib_expected_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin1Distrib.md));

- the third and fourth derivatives in \\\theta\\, whose own cancellation
  in the powers of \\r\\ is untouched
  ([`distrib_deriv3.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin1Distrib.md),
  [`distrib_deriv4.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBin1Distrib.md)).

A fit reaches that regime routinely: on 2,000 Poisson counts with mean 4
the dispersion is estimated at about \\1.7\times 10^{-8}\\ and the run
reports convergence.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Neither estimate is
closed form. The method of moments supplies the starting values
\\\hat\mu = \bar y\\ and \\\hat\theta = s^2/\bar y - 1\\, with \\s^2\\
the sample variance.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \> 0\\ the mean and
\\\theta \> 0\\ the dispersion. \\r = \mu/\theta\\ is the size, \\\psi\\
the digamma function and \\\psi'\\ the trigamma. \\\eta\\ is a parameter
on the unconstrained scale of its link, with \\\theta_j =
g^{-1}(\eta_j)\\.

## References

Cameron, A. C. and Trivedi, P. K. (1986). Econometric models based on
count data: comparisons and applications of some estimators and tests.
*Journal of Applied Econometrics* **1**, 29-53.

## See also

[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the quadratic-variance family;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
for the limit as \\\theta\\ goes to zero and
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md)
for the case \\\mu = \theta\\;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
for a Poisson mixed over an inverse Gaussian;
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for counts with excess zeros;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[NegBin1Distrib](https://statmodels7.github.io/distributions7/reference/NegBin1Distrib.md)
for the class.

## Examples

``` r
d <- negbin1_distrib()
d
#> Distribution: Negbin1
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (mean)               | Link: log        | Domain: (0, Inf)
#>   theta (dispersion)         | Link: log        | Domain: (0, Inf)

# The two negative binomials are different families: at the same (mu, theta)
# this one has variance mu(1 + theta) = 20 and the other mu + mu^2/theta = 8.
th <- list(mu = 4, theta = 4)
c(nb1 = variance(d, th), nb2 = variance(negbin2_distrib(), th))
#> nb1 nb2 
#>  20   8 

# The variance-to-mean ratio is 1 + theta at every mean.
vapply(c(1, 10, 100),
       function(m) variance(d, list(mu = m, theta = 4)) / m, numeric(1))
#> [1] 5 5 5

# At mu = theta the size is 1 and the law is the geometric.
all.equal(distrib_pdf(d, 0:4, th), dgeom(0:4, prob = 1 / 5))
#> [1] TRUE

# As theta goes to zero it is the Poisson.
rbind(nb1 = distrib_pdf(d, 0:4, list(mu = 3, theta = 1e-6)),
      poisson = dpois(0:4, 3))
#>               [,1]      [,2]      [,3]      [,4]      [,5]
#> nb1     0.04978714 0.1493613 0.2240418 0.2240417 0.1680313
#> poisson 0.04978707 0.1493612 0.2240418 0.2240418 0.1680314

# Fitting recovers the parameters; the moment estimates start it off.
set.seed(5)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
rbind(fitted  = coef(fit),
      moments = c(mu = mean(z), theta = var(z) / mean(z) - 1))
#>            mu    theta
#> fitted  3.968 3.991133
#> moments 3.968 3.909959

# On an equidispersed sample the dispersion runs to its boundary, which is
# the regime the expected information page warns about.
set.seed(4)
coef(fit_distrib(d, rpois(2000, 4)))
#>           mu        theta 
#> 3.941000e+00 1.666388e-08 
```
