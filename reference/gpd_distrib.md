# Generalized Pareto Distribution Object

Builds a generalized Pareto distribution object with scale \\\sigma \>
0\\ and shape \\\xi\\. It is the family of exceedances over a high
threshold, and the natural companion of
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
in an analysis of extremes.

At \\\xi = 0\\ it is the exponential, at \\\xi \> 0\\ it has a
polynomial tail, and at \\\xi \< 0\\ it has a **finite upper endpoint**
at \\-\sigma/\xi\\. That last case makes the family non-regular, and it
is the one thing to know before using it.

## Usage

``` r
gpd_distrib(link_sigma = log_link(), link_xi = identity_link())
```

## Arguments

- link_sigma:

  A `linkfunctions7` link object for the scale \\\sigma\\, which must be
  strictly positive. Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_xi:

  A link object for the shape \\\xi\\, which is unconstrained. Defaults
  to
  [`linkfunctions7::identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.html),
  the shape being free to take either sign.

## Value

An S7 object of class
[GPDDistrib](https://statmodels7.github.io/distributions7/reference/GPDDistrib.md),
inheriting from `continuous_distrib`. Its `params` are
`c("sigma", "xi")`, its `bounds` `c(0, Inf)`, and its `link_params` the
two links given here.

## Density and distribution function

\$\$f(y) = \dfrac{1}{\sigma}\left(1 + \dfrac{\xi
y}{\sigma}\right)^{-1/\xi-1}, \qquad F(q) = 1 - \left(1 + \dfrac{\xi
q}{\sigma}\right)^{-1/\xi}.\$\$ At \\\xi = 0\\ both reduce to the
exponential. The implementation reaches that limit through a series, so
the parameter may pass through zero during a fit. The quantile function
is elementary, \\Q(p) = \sigma\\(1-p)^{-\xi}-1\\/\xi\\, so the generator
is an exact inverse transform.

## It is not parametrized by its mean

The mean is \\\sigma/(1-\xi)\\ and exists only for \\\xi \< 1\\, so a
mean parametrization would leave the family undescribable exactly where
it is most used, in the heavy-tailed regime. The same argument keeps
[`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.md)
and
[`variance()`](https://statmodels7.github.io/distributions7/reference/variance.md)
apart for the multivariate \\t\\: a parametrization must not depend on a
moment that need not exist.

## The support depends on the parameters

For \\\xi \< 0\\ the support is \\\[0, -\sigma/\xi\]\\, and this is the
first family here of which that is true. What it costs is the automatic
license to differentiate under the integral sign, on which the Bartlett
identities rest. Two things survive and one does not:

- the derivatives returned are correct as derivatives of the log-density
  at every admissible point, whatever the sign of \\\xi\\;

- the expected information exists and is Smith's closed form for \\\xi
  \> -1/2\\. The condition is exactly that the integrand be integrable:
  near the upper endpoint the second derivative grows like
  \\(1-u)^{-2\|\xi\|}\\ on the probability scale, which is integrable if
  and only if \\\|\xi\| \< 1/2\\;

- below \\\xi = -1/2\\ the information does not exist,
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
  returns `NA`, and the classical asymptotics of the maximum likelihood
  estimator do not hold (Smith, 1985).

The object's `bounds` are `c(0, Inf)` because they are fixed at
construction while the true endpoint moves with the parameters. The
density is zero beyond it, so nothing computes a wrong number, and a
caller reading `bounds` should ask
[`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md)
instead.

## Moments

The mean is \\\sigma/(1-\xi)\\ for \\\xi \< 1\\, the variance
\\\sigma^2/\\(1-\xi)^2(1-2\xi)\\\\ for \\\xi \< 1/2\\, the skewness
needs \\\xi \< 1/3\\ and the kurtosis \\\xi \< 1/4. \\ Each returns
`Inf` above its threshold, so a fitted object can report a mean where it
has no variance.

## Parameter domains

- \\\sigma \in (0, \infty)\\

- \\\xi \in (-\infty, \infty)\\

## References

Smith, R. L. (1985). Maximum likelihood estimation in a class of
nonregular cases. *Biometrika* 72, 67-90.

Davison, A. C. and Smith, R. L. (1990). Models for exceedances over high
thresholds. *Journal of the Royal Statistical Society B* 52, 393-442.

## See also

[`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)
for the \\\xi = 0\\ case,
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md)
and
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
for the other extreme-value families,
[`gpd_endpoint()`](https://statmodels7.github.io/distributions7/reference/gpd_endpoint.md)
for the moving support, and
[GPDDistrib](https://statmodels7.github.io/distributions7/reference/GPDDistrib.md)
for the class and its method list.

## Examples

``` r
d <- gpd_distrib()
d@params
#> [1] "sigma" "xi"   
th <- list(sigma = 1.5, xi = 0.3)

distrib_pdf(d, c(0.2, 1, 4), th)
#> [1] 0.5624677 0.3025450 0.0522069

# Shape zero is the exponential, reached by a series.
all.equal(distrib_pdf(d, c(0.2, 1, 4), list(sigma = 1.5, xi = 0)),
          dexp(c(0.2, 1, 4), rate = 1 / 1.5))
#> [1] TRUE

# Which moments exist depends on the shape, and the mean can be finite
# where the variance is not.
t(vapply(c(-0.3, 0.2, 0.4, 0.6), function(x) {
  p <- list(sigma = 1.5, xi = x)
  c(xi = x, mean = mean(d, p), variance = variance(d, p))
}, numeric(3)))
#>        xi     mean   variance
#> [1,] -0.3 1.153846  0.8321006
#> [2,]  0.2 1.875000  5.8593750
#> [3,]  0.4 2.500000 31.2500000
#> [4,]  0.6 3.750000        Inf

# The information exists only above -1/2.
distrib_expected_hessian(d, 0, list(sigma = 1.5, xi = -0.7))
#> $sigma_sigma
#> [1] NA
#> 
#> $sigma_xi
#> [1] NA
#> 
#> $xi_xi
#> [1] NA
#> 

# A fit recovers both parameters.
set.seed(44)
x <- distrib_rng(d, 4000, th)
coef(fit_distrib(d, x))
#>     sigma        xi 
#> 1.4450802 0.2969134 
```
