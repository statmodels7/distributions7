# Zero-Inflated Distribution Object

Wraps a DISCRETE distribution into a mixture that keeps the parent
intact and adds a second source of zeros, with probability \\\zeta\\
carried by the new parameter `zi`: \$\$P(Y = y; \theta, \zeta) =
\begin{cases} \zeta + (1 - \zeta)f(0; \theta) & y = 0 \\ (1 - \zeta)f(y;
\theta) & y \> 0. \end{cases}\$\$ It is the wrapper to reach for when
the data carry MORE zeros than the parent can produce and a zero may
plausibly have come either from the count process or from a mechanism
that switches it off.

## Usage

``` r
zero_inflated(distrib, link_zi = logit_link())
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` whose support includes 0,
  such as
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
  or
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).
  A continuous distribution, an already-wrapped one, and a support of
  fewer than `n_params + 2` points are each rejected with an error
  saying which condition failed.

- link_zi:

  The link carrying \\\zeta\\ to the unconstrained scale, a
  [`linkfunctions7::link`](https://statmodels7.github.io/linkfunctions7/reference/link.html)
  object. Defaults to
  [`linkfunctions7::logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  which keeps it strictly inside \\(0, 1)\\ at every point of the free
  scale.

## Value

An S7 object of class
[ZeroInflatedDistrib](https://statmodels7.github.io/distributions7/reference/ZeroInflatedDistrib.md),
inheriting from `discrete_distrib`. Its `params` are the parent's
followed by `zi`; `n_params` is the parent's plus one;
`params_bounds$zi` is \\(0, 1)\\; `link_params$zi` is `link_zi` and the
rest are the parent's; `params_interpretation` gains
`"prob. of structural zero"`; and `distrib_name` is `"zero-inflated "`
followed by the parent's.

## Zero-inflation against zero-adjustment

The two wrappers differ in what they do to the mass the parent already
places at zero. Zero-inflation ADDS to it, so \\P(Y = 0) = \zeta +
(1-\zeta)f(0) \> f(0)\\: the model can produce more zeros than the
parent and never fewer, and the observed zeros are a mixture of
structural and sampling ones that no single observation can be assigned
to. Zero-adjustment REPLACES it: the parent is truncated away from zero
and the mass there becomes a free parameter, which may sit above or
below \\f(0)\\.

A hurdle model therefore also handles UNDER-dispersed zeros, and its
likelihood factorizes into a binary part and a positive-count part that
can be read separately. Zero-inflation keeps the parent's
interpretation, so \\\theta\\ still describes the count process the
non-structural observations come from, while the hurdle re-interprets
\\\theta\\ as the parameters of a truncated law.

## What the parent must be

Zero-inflation adds mass to a zero that already carries some, so it
needs a discrete distribution with 0 in its support. A continuous one
has \\P(Y = 0) = 0\\ and nothing to inflate; placing a point mass at
zero beside a density is zero-ADJUSTMENT, and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
handles it. The constructor also fails where the result would not be
identified:

- the parent already models a probability of zero, since only the total
  mass at zero would then be identified;

- the support is too small for one more parameter. A distribution on
  \\k\\ points has \\k-1\\ free probabilities, so at least
  `n_params + 2` support points are needed. This rules out the Bernoulli
  and `binomial_distrib(size = 1)`, where `mu` and `zi` between them
  describe a single free cell.

A large support is NECESSARY without being sufficient. With \\\mu\\
large enough that \\f(0)\\ underflows, or \\\zeta\\ close to zero, the
ridge reappears in the data instead of in the model, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reports the standard errors that reveal it.

## What the result supports

The whole `distrib` contract: the mass function, the distribution
function, the quantile function, the generator, the analytic score, the
observed and the expected Hessian, and the third and fourth derivatives
from the shared wrapper machinery. The moments come from
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md)
numerically.

## Notation

\\f\\ is the parent's mass function, \\\theta\\ its parameters,
\\\zeta\\ the probability of a structural zero and \\k\\ the number of
support points.

## See also

[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the hurdle counterpart,
[ZeroInflatedDistrib](https://statmodels7.github.io/distributions7/reference/ZeroInflatedDistrib.md)
for the class,
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
which the hurdle uses internally, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate the result.

## Examples

``` r
zip <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.2)
zip@params
#> [1] "mu" "zi"

distrib_pdf(zip, 0:5, theta)
#> [1] 0.23982965 0.11948896 0.17923345 0.17923345 0.13442508 0.08065505

# More mass at zero than the Poisson alone can put there, and the rest of
# the mass function shrunk by 1 - zi.
c(inflated = distrib_pdf(zip, 0, theta), poisson = dpois(0, 3))
#>   inflated    poisson 
#> 0.23982965 0.04978707 
all.equal(distrib_pdf(zip, 1:5, theta), 0.8 * dpois(1:5, 3))
#> [1] TRUE

# It is still a distribution: the mass sums to one.
sum(distrib_pdf(zip, 0:200, theta))
#> [1] 1

# A fit recovers both parameters, the extra zeros identifying zi.
set.seed(1)
y <- distrib_rng(zip, 2000, theta)
round(coef(fit_distrib(zip, y)), 3)
#>    mu    zi 
#> 2.970 0.211 

# A Bernoulli has no room for a second parameter, and a wrapper cannot be
# stacked on another.
try(zero_inflated(bernoulli_distrib()))
#> Error : zero_inflated() cannot wrap 'bernoulli': its support has 2 points, so the family has 1 free
#>   probability, while the wrapped distribution would have 2 parameters. They are not
#>   identified -- different parameter values give exactly the same distribution.
#>   A support of at least 3 points is required.
try(zero_inflated(zero_adjusted(poisson_distrib())))
#> Error : zero_inflated() cannot wrap 'zero-adjusted poisson', which already models the probability of a zero.
#>   Stacking the two leaves only their combination identified: the second
#>   parameter has an identically zero score, and any optimizer will wander
#>   along that ridge. Apply exactly one zero wrapper to a plain distribution.

# Nor can a continuous parent be inflated: there is no mass at zero to add
# to, and the message names the wrapper that does apply.
try(zero_inflated(gaussian1_distrib()))
#> Error : zero_inflated() requires a discrete distribution: a continuous one has
#>   P(Y = 0) = 0, so there is no mass at zero to inflate. To place a point mass
#>   at zero alongside a density, use zero_adjusted().
```
