# Beta-Binomial Distribution, Two Shapes

Builds the distribution object for the beta-binomial family in its
canonical parametrization, the two beta shapes \\\alpha \> 0\\ and
\\\beta \> 0\\, on the finite support \\\\0, 1, \dots, n\\\\. The
returned object carries closed-form derivatives of the log-mass to
fourth order and expectations that are exact finite sums over the
support, so every generic of the toolkit answers without a quadrature.

This is the same law as
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
which is written in a mean proportion and a dispersion. Use this one
when the shapes themselves are the quantities of interest, and that one
when a regression on the mean is.

## Usage

``` r
betabinom2_distrib(size, link_alpha = log_link(), link_beta = log_link())
```

## Arguments

- size:

  The number of trials \\n\\, a single positive integer. It is a
  constant of the distribution and not a parameter, as for
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
  so an object cannot be reused across data sets whose group sizes
  differ. Anything else signals an error naming the argument.

- link_alpha:

  A `link` object from `linkfunctions7` for the shape \\\alpha\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

- link_beta:

  A `link` object from `linkfunctions7` for the shape \\\beta\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  for the same reason.

## Value

An S7 object of class `BetaBinom2Distrib`, inheriting from
`discrete_distrib`, with `size` the trial count, `distrib_name`
`"betabinom2 [size=n]"`, `dimension` `"univariate"`, `bounds`
`c(0, size)`, `params` `c("alpha", "beta")`, `n_params` `2`,
`params_bounds` the domain \\(0, \infty)\\ for both, and `link_params`
the two links given here.

## The parametrization

The mass on \\y \in \\0, \dots, n\\\\ is \$\$P(Y=y) =
\binom{n}{y}\frac{B(y+\alpha,\\ n-y+\beta)}{B(\alpha, \beta)},\$\$ with
\\B\\ the beta function, and the family is the binomial with its success
probability drawn from \\\mathrm{Beta}(\alpha, \beta)\\. Writing \\S =
\alpha + \beta\\ for the concentration, the moments are
\$\$\mathbb{E}\[Y\] = \frac{n\alpha}{S}, \qquad \operatorname{Var}(Y) =
\frac{n\alpha\beta\\(S+n)}{S^{2}(S+1)}.\$\$ The correspondence with
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
is \\\mu = \alpha/S\\ and \\\sigma = 1/S\\, so a large concentration is
a small dispersion and the binomial limit.

## Derivatives

The log-mass is a sum of log-gamma terms, so a derivative of order \\k\\
replaces each by \\\psi^{(k-1)}\\ and all four orders come from one
routine,
[`betabinom2_derivs()`](https://statmodels7.github.io/distributions7/reference/betabinom2_derivs.md).
A component naming both shapes keeps only the terms in \\S\\ and is
therefore free of the data at every order.

Every expectation is an **exact finite sum** over \\\\0, \dots, n\\\\:
the support is bounded, so there is nothing to integrate over.

## The cancellation at a large concentration

The two beta functions of the mass are each of magnitude \\S\log S\\
while their difference is of order one, so writing the mass as that
difference loses one digit per factor of ten in \\S\\.
[`betabinom_log_mass()`](https://statmodels7.github.io/distributions7/reference/betabinom_log_mass.md)
switches to a sum of logarithms past a measured threshold and stays
exact: at \\n = 10\\, \\y = 3\\ and \\\alpha/S = 0.4\\, the log-mass
agrees with the binomial one to twelve figures at \\S = 10^{14}\\, where
the direct route is wrong in the third decimal.

**The derivatives are not rewritten that way here**, and cede earlier
than the mass. Measured at the same setting, the relative error of the
\\\alpha\\ score against its limit is \\6\times10^{-8}\\ at \\S =
10^8\\, \\4\times10^{-4}\\ at \\10^{12}\\ and 1.8 at \\10^{15}\\.
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)'s
compiled kernel forms the same differences as sums of reciprocals and
holds nine figures at \\S = 10^{15}\\, so that is the parametrization to
use where such a concentration is reachable.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Neither shape is closed
form, and the two are strongly correlated at a large concentration, the
data then determining their ratio far better than their size.

## Notation

\\\ell\\ is the log-mass of one observation, \\\alpha, \beta \> 0\\ the
two beta shapes, \\S = \alpha+\beta\\ the concentration, \\n\\ the trial
count, \\B\\ the beta function and \\\psi^{(k)}\\ the polygamma
functions. \\\eta\\ is a parameter on the unconstrained scale of its
link, with \\\theta = g^{-1}(\eta)\\.

## References

Skellam, J. G. (1948). A probability distribution derived from the
binomial distribution by regarding the probability of success as
variable between the sets of trials. *Journal of the Royal Statistical
Society, Series B*, **10**(2), 257-261.

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). *Univariate Discrete
Distributions*, 3rd edition, Section 6.9. Wiley, Hoboken.

## See also

[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
for the same law in a mean proportion and a dispersion, which is the one
to model a mean with;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
for the limit at a large concentration and
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for the mixing law;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the shapes;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[BetaBinom2Distrib](https://statmodels7.github.io/distributions7/reference/BetaBinom2Distrib.md)
for the class.

## Examples

``` r
d <- betabinom2_distrib(size = 10)
d
#> Distribution: Betabinom2 [size=10]
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   alpha (shape)              | Link: log        | Domain: (0, Inf)
#>   beta  (shape)              | Link: log        | Domain: (0, Inf)

# The mass over the support sums to one.
th <- list(alpha = 2, beta = 3)
sum(distrib_pdf(d, 0:10, th))
#> [1] 1

# The closed-form moments, and the same numbers from the shapes.
c(mean = mean(d, th), var = variance(d, th),
  closed_mean = 10 * 2 / 5, closed_var = 10 * 2 * 3 * 15 / (25 * 6))
#>        mean         var closed_mean  closed_var 
#>           4           6           4           6 

# The same law as betabinom1 at mu = alpha / S and sigma = 1 / S.
all.equal(distrib_pdf(d, 0:10, th),
          distrib_pdf(betabinom1_distrib(size = 10), 0:10,
                      list(mu = 0.4, sigma = 0.2)))
#> [1] TRUE

# At a large concentration the mass is still exact against the binomial
# limit, where forming it from two beta functions is not.
S <- 1e14
c(shipped = distrib_pdf(d, 3, list(alpha = 0.4 * S, beta = 0.6 * S),
                        log = TRUE),
  two_betas = lchoose(10, 3) +
    (lbeta(3 + 0.4 * S, 7 + 0.6 * S) - lbeta(0.4 * S, 0.6 * S)),
  binomial = dbinom(3, 10, 0.4, log = TRUE))
#>   shipped two_betas  binomial 
#> -1.537160 -1.540633 -1.537160 

# Fitting recovers both shapes.
set.seed(5)
z <- distrib_rng(d, 4000, th)
coef(fit_distrib(d, z))
#>    alpha     beta 
#> 2.005073 2.992054 
```
