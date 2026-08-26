# Beta-Binomial Distribution, Mean Proportion and Dispersion

Builds the distribution object for the beta-binomial family parametrized
by a mean proportion \\\mu \in (0, 1)\\ and a dispersion \\\sigma \>
0\\, on the finite support \\\\0, 1, \dots, n\\\\. The returned object
carries closed-form derivatives of the log-mass to fourth order and
expectations that are exact finite sums over the support, so every
generic of the toolkit answers without a quadrature.

The family is the binomial with its success probability drawn from a
beta, which makes it the standard model for a proportion whose trials
are not independent. It is overdispersed relative to a binomial at every
\\\sigma \> 0\\ and approaches one as \\\sigma \to 0\\.

## Usage

``` r
betabinom1_distrib(size, link_mu = logit_link(), link_sigma = log_link())
```

## Arguments

- size:

  The number of trials \\n\\, a single positive integer. It is a
  constant of the distribution and not a parameter, as for
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
  so an object cannot be reused across data sets whose group sizes
  differ. Anything else signals an error naming the argument.

- link_mu:

  A `link` object from `linkfunctions7` for the mean proportion \\\mu\\.
  Defaults to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  which maps \\(0, 1)\\ onto the line.

- link_sigma:

  A `link` object from `linkfunctions7` for the dispersion \\\sigma\\.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html),
  which maps \\(0, \infty)\\ onto the line and so keeps every fitted
  value positive.

## Value

An S7 object of class `BetaBinom1Distrib`, inheriting from
`discrete_distrib`, with `size` the trial count, `distrib_name`
`"beta-binomial [size=n]"`, `dimension` `"univariate"`, `bounds`
`c(0, size)`, `params` `c("mu", "sigma")`, `n_params` `2`,
`params_bounds` the domains \\(0, 1)\\ and \\(0, \infty)\\, and
`link_params` the two links given here.

## The parametrization

Write \\\alpha = \mu/\sigma\\ and \\\beta = (1-\mu)/\sigma\\ for the two
beta shapes. The mass on \\y \in \\0, \dots, n\\\\ is \$\$P(Y = y) =
\binom{n}{y} \dfrac{B(y+\alpha,\\ n-y+\beta)}{B(\alpha, \beta)},\$\$
with \\B\\ the beta function.
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
is the same law written in \\(\alpha, \beta)\\ directly; the two agree
mass for mass at \\\mu = \alpha/(\alpha+\beta)\\ and \\\sigma =
1/(\alpha+\beta)\\.

The family is **not** reachable from anything else in the package. It is
neither a binomial with a parameter held fixed nor a wrapper over one,
the mixing being over the success probability and not over the outcome.

## Overdispersion

The mean is \\n\mu\\ and the variance \$\$\operatorname{Var}(Y) =
n\mu(1-\mu) \left(1 + (n-1)\dfrac{\sigma}{1+\sigma}\right),\$\$ so the
inflation over a binomial of the same mean is \\1 +
(n-1)\sigma/(1+\sigma)\\, always above 1 and rising with both the trial
count and the dispersion. The factor \\\sigma/(1+\sigma)\\ is the
**intraclass correlation** between two trials of the same group. At \\n
= 10\\ and \\\sigma = 0.5\\ that correlation is \\1/3\\ and the variance
is four times a binomial's, 8.4 against 2.1.

As \\\sigma \to 0\\ the mass converges to the binomial's at rate
\\O(\sigma)\\, and the compiled kernel stays accurate all the way there;
see the cancellation note below.

## Derivatives

The parameters enter only through the two shapes, where every derivative
of order \\k\\ is a difference of \\\psi^{(k-1)}\\: \$\$\dfrac{\partial
\ell}{\partial \alpha} = \psi(y+\alpha) - \psi(\alpha) - \psi(n+S) +
\psi(S), \qquad S = \alpha + \beta,\$\$ and likewise in \\\beta\\. The
reported components follow by the chain rule of \\(\alpha, \beta) =
(\mu/\sigma, (1-\mu)/\sigma)\\, which is linear in \\\mu\\ at fixed
\\\sigma\\, so every partial of the map carrying two or more \\\mu\\ is
exactly zero and the third and fourth orders are a short partition sum.
Orders three and four are therefore **closed form** as well, in
[`distrib_deriv3.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom1Distrib.md),
and not the numerical fallback.

Every expectation is an **exact finite sum** over \\\\0, \dots, n\\\\:
the support is bounded, so there is nothing to integrate over.

## The cancellation at a small dispersion

The two beta functions of the mass are of magnitude
\\(\alpha+\beta)\log(\alpha+\beta)\\ while their difference is of order
one, so writing the mass as that difference loses one digit per factor
of ten in the concentration \\1/\sigma\\. At \\\sigma = 10^{-14}\\ the
direct route is wrong in the third decimal of the log-mass. The kernel
switches to a sum of logarithms instead, the shifts being integers, and
the same rewrite runs through the score: measured at \\n = 10\\, \\\mu =
0.4\\, \\y = 3\\, the mean component holds the binomial value \\-25/6\\
to nine figures at a concentration of \\10^{15}\\.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood on the link scale. Neither estimate is
closed form. A sample with no extra-binomial variation drives
\\\hat\sigma\\ towards zero, which is the correct answer and the
boundary of the parameter space.

## Notation

\\\ell\\ is the log-mass of one observation, \\\mu \in (0,1)\\ the mean
proportion, \\\sigma \> 0\\ the dispersion, \\n\\ the trial count, \\B\\
the beta function and \\\psi\\ the digamma function. \\\eta\\ is a
parameter on the unconstrained scale of its link, with \\\theta =
g^{-1}(\eta)\\.

## References

Skellam, J. G. (1948). A probability distribution derived from the
binomial distribution by regarding the probability of success as
variable between the sets of trials. *Journal of the Royal Statistical
Society, Series B*, **10**(2), 257-261.

Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). *Univariate Discrete
Distributions*, 3rd edition, Section 6.9. Wiley, Hoboken.

## See also

[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for the same law in its two beta shapes;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
for the limit at \\\sigma \to 0\\ and
[`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)
for the mixing law;
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the unbounded-count analogue, a Poisson mixed over a gamma;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the parameters;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[BetaBinom1Distrib](https://statmodels7.github.io/distributions7/reference/BetaBinom1Distrib.md)
for the class.

## Examples

``` r
d <- betabinom1_distrib(size = 10)
d
#> Distribution: Beta-binomial [size=10]
#> Type:         Discrete
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu    (mean proportion)    | Link: logit      | Domain: (0, 1)
#>   sigma (dispersion)         | Link: log        | Domain: (0, Inf)

# The mass over the support sums to one.
th <- list(mu = 0.3, sigma = 0.5)
sum(distrib_pdf(d, 0:10, th))
#> [1] 1

# Four times the variance of a binomial with the same mean: the intraclass
# correlation is sigma / (1 + sigma) = 1/3, and 1 + 9/3 = 4.
c(mean = mean(d, th), var = variance(d, th),
  binomial_var = variance(binomial_distrib(size = 10), list(mu = 0.3)),
  icc = 0.5 / 1.5)
#>         mean          var binomial_var          icc 
#>    3.0000000    8.4000000    2.1000000    0.3333333 

# The dispersion going to zero is the binomial, at rate O(sigma).
vapply(c(1e-2, 1e-4, 1e-6), function(s)
  max(abs(distrib_pdf(d, 0:10, list(mu = 0.3, sigma = s)) -
          dbinom(0:10, 10, 0.3))), numeric(1))
#> [1] 1.248997e-02 1.333223e-04 1.334134e-06

# The same law as betabinom2 at the implied shapes.
all.equal(distrib_pdf(d, 0:10, list(mu = 0.4, sigma = 0.2)),
          distrib_pdf(betabinom2_distrib(size = 10), 0:10,
                      list(alpha = 2, beta = 3)))
#> [1] TRUE

# Fitting recovers both parameters.
set.seed(3)
z <- distrib_rng(d, 2000, list(mu = 0.35, sigma = 0.4))
coef(fit_distrib(d, z))
#>        mu     sigma 
#> 0.3568109 0.4207828 
```
