# Multinomial Distribution

Builds the distribution object for the multinomial family on \\p\\
categories and \\n\\ trials, parametrized by a probability vector
carried on a `parameters7` simplex. The returned object carries
closed-form derivatives of the log-mass to fourth order and a
closed-form expected information, and its support is small enough to
enumerate, so every expectation can be checked as an exact sum.

## Usage

``` r
multinomial_distrib(n_dim, size, probs = parameters7::simplex(n_dim))
```

## Arguments

- n_dim:

  The number of categories \\p\\, a single integer of at least 2.
  Anything else signals an error naming the argument.

- size:

  The number of trials \\n\\, a single positive integer. It is a
  constant of the distribution and not a parameter, as for
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
  so an object cannot be reused across data sets whose trial counts
  differ.

- probs:

  A `parameters7` parameter producing \\p\\ coordinates that sum to one,
  normally a
  [`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html).
  Defaults to `parameters7::simplex(n_dim)`. An object that is not a
  `parameters7` parameter, or that produces a different number of
  coordinates, signals an error. Its free values are unconstrained
  already and carry the identity link.

## Value

An S7 object of class `MultinomialDistrib`, inheriting from
`multivariate_distrib`, with `size` the trial count, `param` the simplex
given here, `distrib_name` `"multinomial [pd, size=n, probs=<chart>]"`,
`dimension` `"multivariate"`, `n_dim` \\p\\, `bounds` `c(0, size)`,
`params` the simplex's free names prefixed by `probs_`, `n_params`
\\p-1\\, and `link_params` the identity for each.

## The parametrization

The mass on the weak compositions of \\n\\ into \\p\\ parts is \$\$P(Y =
y) = \frac{n!}{\prod\_{j=1}^{p} y_j!}\prod\_{j=1}^{p} p_j^{y_j},\$\$
with \$\$\mathbb{E}\[Y_j\] = n p_j, \qquad \operatorname{Var}(Y_j) = n
p_j (1 - p_j), \qquad \operatorname{Cov}(Y_j, Y_k) = -n p_j p_k.\$\$
There is **no dispersion parameter**: the probabilities are all the
family has, and \\p-1\\ of them are free.

The probabilities are carried by a `parameters7` simplex and flattened
into scalars with identity links, exactly as a covariance is for the
multivariate gaussian. The constraint that they be positive and sum to
one lives in the parameter, where a scalar link could not express it.

## A multivariate family that is discrete

This is the first family of the package that is multivariate and
discrete, so its support is a finite set of points rather than a region.
[`mv_support.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_support.MultinomialDistrib.md)
enumerates them, which lets an expectation be an exact sum and lets
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
test the total mass by addition: measured, the mass over the support
comes back at 1 to `8e-16` and the closed-form information agrees with
the summed observed Hessian to the same order. An importance-sampling
check would only ever compare against Monte Carlo error, and a
normalization wrong by a thousandth would pass it.

The count of support points is \\\binom{n+p-1}{p-1}\\ and grows quickly,
so the enumeration is a validation tool rather than a fitting route.

## Score and information

With \\A = \partial p/\partial\eta\\ the chart's Jacobian,
\$\$\dfrac{\partial\ell}{\partial\eta_k} = \sum_j
\dfrac{y_j}{p_j}A\_{jk}, \qquad \mathbb{E}\[\ell^{(\eta_k\eta_l)}\] =
-n\sum_j \dfrac{A\_{jk}A\_{jl}}{p_j}.\$\$ The expected form is closed
because \\\mathbb{E}\[y_j\] = np_j\\ turns the second-derivative term
into \\n\sum_j B\_{j,kl}\\, which vanishes: the probabilities sum to
one, so every derivative of their sum is zero.

On the default additive log-ratio chart the resulting information is
exactly the covariance of the first \\p-1\\ counts, since \\A\_{jk} =
p_j(\delta\_{jk}-p_k)\\ collapses the sum to \\\delta\_{kl}p_k -
p_kp_l\\. Measured, the two agree to `2e-16` at every dimension and
trial count tried.

## Marginals and conjugacy

Coordinate \\j\\ is \\\mathrm{Binomial}(n, p_j)\\, the other categories
collapsing into a single failure, so
[`mv_marginal.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MultinomialDistrib.md)
returns an object. The Dirichlet is the conjugate prior for the
probability vector, and
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
is written on the same simplex.

## Estimation

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the log-likelihood over the simplex's free values, whose links
are the identity, so the two scales coincide. The maximum likelihood
estimate of \\p_j\\ is the pooled sample proportion.

## Notation

\\\ell\\ is the log-mass of one observation, \\p\\ the probability
vector, \\n\\ the trial count, \\\eta\\ the free vector of the simplex
chart, \\A = \partial p/\partial\eta\\ its Jacobian and \\B\\ its
second-derivative arrays.

## References

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1997). *Discrete
Multivariate Distributions*, Chapter 35. Wiley, New York.

Agresti, A. (2013). *Categorical Data Analysis*, 3rd edition, Section
1.2. Wiley, Hoboken.

## See also

[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
for a category's marginal and the two-category case;
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
for the conjugate family on the same simplex;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
for the overdispersed two-category count;
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
for the chart the probabilities ride;
[`numericals7::compositions()`](https://statmodels7.github.io/numericals7/reference/compositions.html)
for the support enumeration;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
to estimate the probabilities;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate a family of your own against the same battery this one
passes;
[MultinomialDistrib](https://statmodels7.github.io/distributions7/reference/MultinomialDistrib.md)
for the class.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
d
#> Distribution: Multinomial [3d, Size=5, Probs=simplex]
#> Type:         Continuous, 3-dimensional
#> Dimensions:   multivariate
#> 
#> Parameters:
#>   probs_alr1 (probability)        | Link: identity   | Domain: (-Inf, Inf)
#>   probs_alr2 (probability)        | Link: identity   | Domain: (-Inf, Inf)

# Two free probabilities for three categories; no dispersion parameter.
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
mv_location(d, th)
#> [1] 2.130063 1.291948 1.577989

# The support is a finite set of points, so the mass sums exactly.
supp <- mv_support(d, th)
c(points = nrow(supp), mass = sum(distrib_pdf(d, supp, th)))
#> points   mass 
#>     21      1 

# The covariance is singular, its diagonal the binomial variances.
pr <- mv_location(d, th) / 5
S <- mv_sigma(d, th)
c(rank = qr(S)$rank, dim = 3)
#> rank  dim 
#>    2    3 
rbind(diagonal = diag(S), binomial = 5 * pr * (1 - pr))
#>              [,1]      [,2]     [,3]
#> diagonal 1.222629 0.9581222 1.079979
#> binomial 1.222629 0.9581222 1.079979

# On the default chart the expected information is the covariance of the
# first p - 1 counts.
eh <- distrib_expected_hessian(d, matrix(0, 1, 3), th)
rbind(information = -vapply(eh, function(v) v[1], numeric(1)),
      covariance = c(S[1, 1], S[2, 2], S[1, 2]))
#>             probs_alr1_probs_alr1 probs_alr2_probs_alr2 probs_alr1_probs_alr2
#> information              1.222629             0.9581222            -0.5503861
#> covariance               1.222629             0.9581222            -0.5503861

# Fitting recovers the probabilities.
set.seed(4)
coef(fit_distrib(d, distrib_rng(d, 1500, th)))
#> probs_alr1 probs_alr2 
#>  0.2708608 -0.2027115 
```
