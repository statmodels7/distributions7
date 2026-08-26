# Beta Distribution, the Two Shapes

Builds the distribution object for the beta family on \\(0, 1)\\ in its
canonical parametrization, the two shapes \\\alpha \> 0\\ and \\\beta \>
0\\. The returned object carries closed-form derivatives of the
log-density to fourth order, in the parameters and in the response, and
closed-form moments, so every generic of the toolkit answers without a
numerical fallback.

The two arguments choose the links that carry each parameter to the
unconstrained scale an optimizer works on. Both default to the
logarithm, both shapes being positive.

## Usage

``` r
beta2_distrib(link_alpha = log_link(), link_beta = log_link())
```

## Arguments

- link_alpha:

  A `link` object from `linkfunctions7` for the first shape \\\alpha\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

- link_beta:

  A `link` object from `linkfunctions7` for the second shape \\\beta\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `Beta2Distrib`, inheriting from
`continuous_distrib`, with `distrib_name` `"beta2"`, `dimension`
`"univariate"`, `bounds` `c(0, 1)`, `params` `c("alpha", "beta")`,
`n_params` `2`, `params_bounds` the domain \\(0, \infty)\\ for both, and
`link_params` the two links given here.

## The parametrization

The density on \\y \in (0, 1)\\ is \$\$f(y; \alpha, \beta) =
\dfrac{y^{\alpha-1}(1-y)^{\beta-1}} {B(\alpha, \beta)},\$\$ with \\B\\
the beta function, the distribution function \\F(q) = I_q(\alpha,
\beta)\\ the regularized incomplete beta function, and the quantile
function its numerical inverse. The mean is \\\alpha/(\alpha+\beta)\\
and the variance \\\alpha\beta/\\(\alpha+\beta)^2(\alpha+\beta+1)\\\\.
At \\\alpha = \beta = 1\\ the density is the uniform, and where either
shape falls below one the density is unbounded at the corresponding
endpoint.

This is the same law as
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md),
which carries the mean and a precision, with \\\alpha = \mu\varphi\\ and
\\\beta = (1-\mu)\varphi\\. The mean parametrization is the one a
regression wants; this one is how the family is usually written, and it
is what a conjugate analysis produces, the beta being conjugate for a
binomial probability.

## Where the response stops entering

The beta is an exponential family in the shapes, with sufficient
statistics \\\log y\\ and \\\log(1-y)\\, so the score is
\$\$\dfrac{\partial \ell}{\partial \alpha} = \log y - \psi(\alpha) +
\psi(\alpha+\beta), \qquad \dfrac{\partial \ell}{\partial \beta} =
\log(1-y) - \psi(\beta) + \psi(\alpha+\beta),\$\$ each a statistic minus
its expectation. The data enter only through \\(\alpha-1)\log y +
(\beta-1)\log(1-y)\\, which is linear in the two parameters, so the
second derivative already kills it: \$\$\ell^{(\alpha\alpha)} =
\psi'(\alpha+\beta) - \psi'(\alpha), \qquad \ell^{(\alpha\beta)} =
\psi'(\alpha+\beta), \qquad \ell^{(\beta\beta)} = \psi'(\alpha+\beta) -
\psi'(\beta).\$\$

Three consequences follow. The observed and expected Hessians are the
same matrix, so Fisher scoring and Newton's method take the same step on
the parameter scale, and asking any method for `expected = TRUE` returns
the same numbers. Every third and fourth derivative is likewise a
difference of polygamma functions and free of the data. And the mixed
entry \\\psi'(\alpha+\beta)\\ is positive at every setting, so the two
shapes are never orthogonal.

The derivatives of the *distribution* function in the parameters have no
elementary form, the derivative of an incomplete beta in its shapes
being hypergeometric, and are taken by finite difference on the analytic
cdf.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Neither estimate is
closed form, both estimating equations involving digamma functions. The
method of moments supplies the starting values through the concentration
\\k = \bar y(1-\bar y)/s^2 - 1\\, with \\\hat\alpha = k\bar y\\ and
\\\hat\beta = k(1-\bar y)\\, and the example below shows them landing
beside the estimates.

## Notation

\\\ell\\ is the log-density of one observation and \\\alpha, \beta \>
0\\ the two shapes. \\\psi\\ is the digamma function and \\\psi^{(m)}\\
its \\m\\th derivative. \\\eta\\ is a parameter on the unconstrained
scale of its link, with \\\theta = g^{-1}(\eta)\\.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995). *Continuous
Univariate Distributions*, Volume 2, 2nd edition, Chapter 25. Wiley, New
York.

## See also

[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for the same law in the mean and a precision, which is the
parametrization a regression uses;
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for the beta as a mixing law over a binomial probability;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
for the family a ratio of gammas turns into this one;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[Beta2Distrib](https://statmodels7.github.io/distributions7/reference/Beta2Distrib.md)
for the class.

## Examples

``` r
d <- beta2_distrib()
d
#> Distribution: Beta2
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   alpha (shape)              | Link: log        | Domain: (0, Inf)
#>   beta  (shape)              | Link: log        | Domain: (0, Inf)

# The density is stats::dbeta at these two shapes.
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)
all.equal(distrib_pdf(d, y, th), dbeta(y, 2, 5))
#> [1] TRUE

# Moments: the mean is the ratio, and neither parameter is a mean.
c(mean = mean(d, th), var = variance(d, th),
  skew = skewness(d, th), kurt = kurtosis(d, th))
#>       mean        var       skew       kurt 
#>  0.2857143  0.0255102  0.5962848 -0.1200000 
c(2 / 7, 2 * 5 / (7^2 * 8))
#> [1] 0.2857143 0.0255102

# The same law as beta1 at mu = alpha/(alpha + beta), phi = alpha + beta.
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(beta1_distrib(), y, list(mu = 2 / 7, phi = 7)))
#> [1] TRUE

# Fitting recovers the shapes; the moment estimates start it off.
set.seed(8)
z <- distrib_rng(d, 2000, th)
fit <- fit_distrib(d, z)
m <- mean(z)
k <- m * (1 - m) / var(z) - 1
rbind(fitted  = coef(fit),
      moments = c(alpha = k * m, beta = k * (1 - m)))
#>            alpha     beta
#> fitted  1.927355 4.827894
#> moments 1.933137 4.847228

# Every derivative past the first is free of the response, so Newton and
# Fisher scoring invert the same matrix and reach the same point.
rbind(newton = coef(fit_distrib(d, z, method = "newton")),
      fisher = coef(fit_distrib(d, z, method = "fisher")))
#>           alpha     beta
#> newton 1.927355 4.827894
#> fisher 1.927355 4.827894
```
