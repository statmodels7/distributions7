# Zero-Inflated Distribution Object (Discrete)

Creates a zero-inflated version of an existing **discrete**
distribution: a mixture that keeps the parent intact and adds a second
source of zeros, with probability \\\zeta\\ (parameter `zi`).

Zero-inflation is the right wrapper when the data contain *more* zeros
than the parent can produce, and a zero can plausibly have come either
from the count process or from a separate mechanism that switches it
off. If instead the zeros come from one identifiable mechanism and the
positive values from another, the appropriate model is the hurdle,
[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

## Usage

``` r
zero_inflated(distrib, link_zi = logit_link())
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` whose support includes 0,
  e.g.
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md)
  or
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- link_zi:

  A link function object for the zero-inflation probability \\\zeta\\.
  Defaults to
  [`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html).

## Value

An S7 object of class `ZeroInflatedDistrib` (inheriting from
`discrete_distrib`).

## Details

\$\$ P(Y=y; \theta, \zeta) = \begin{cases} \zeta + (1 - \zeta)f(0;
\theta) & y = 0 \\ (1 - \zeta)f(y; \theta) & y \> 0 \end{cases} \$\$

**Zero-inflation versus zero-adjustment.** The two wrappers differ in
what they do to the mass the parent already places at zero.
Zero-inflation *adds* to it, so \\P(Y = 0) = \zeta + (1-\zeta)f(0) \>
f(0)\\: the model can only ever produce more zeros than the parent,
never fewer, and the observed zeros are a mixture of structural and
sampling ones that no single observation can be assigned to.
Zero-adjustment
([`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md))
*replaces* it: the parent is truncated away from zero and the mass at
zero becomes a free parameter, which can be above or below \\f(0)\\. A
hurdle model therefore also handles *under*-dispersed zeros, and its
likelihood factorises into a binary part and a positive-count part that
can be read separately. Zero-inflation keeps the parent's interpretation
— \\\theta\\ still describes the count process the non-structural
observations come from — while the hurdle re-interprets \\\theta\\ as
the parameters of a truncated law.

**What the parent must be.** Zero-inflation adds mass to a zero that
already carries some, so it requires a discrete distribution with \\0\\
in its support. A continuous distribution has \\P(Y = 0) = 0\\ and
nothing to inflate; putting a point mass at zero next to a density is
zero-*adjustment*, and
[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
handles it. Constructing the object also fails when the result would not
be identified:

- the parent already models a probability of zero (a wrapper cannot be
  stacked on another wrapper: only the total mass at zero would be
  identified);

- the support is too small for one more parameter — a distribution on
  \\k\\ points has \\k-1\\ free probabilities, so at least
  `n_params + 2` support points are needed. This rules out the Bernoulli
  and `binomial_distrib(size = 1)`, where `mu` and `zi` between them
  describe a single free cell.

A large support is necessary but not sufficient: with \\\mu\\ large
enough that \\f(0)\\ underflows, or \\\zeta\\ close to 0, the ridge
reappears in the data rather than in the model.
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reports the standard errors that reveal it.

The resulting object supports the full `distrib` API: pdf, cdf,
quantile, rng, analytical gradient, observed and expected Hessian (all
derived from the parent's), plus numerical moments via
[`moment`](https://statmodels7.github.io/distributions7/reference/moment.md).

## See also

[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the hurdle counterpart,
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
to validate the result.

## Examples

``` r
zip <- zero_inflated(poisson_distrib())
distrib_pdf(zip, 0:5, list(mu = 3, zi = 0.2))
#> [1] 0.23982965 0.11948896 0.17923345 0.17923345 0.13442508 0.08065505

# More mass at zero than the Poisson alone can put there
distrib_pdf(zip, 0, list(mu = 3, zi = 0.2)) > dpois(0, 3)
#> [1] TRUE

# A Bernoulli has no room for a second parameter
try(zero_inflated(bernoulli_distrib()))
#> Error : zero_inflated() cannot wrap 'bernoulli': its support has 2 points, so the family has 1 free
#>   probability, while the wrapped distribution would have 2 parameters. They are not
#>   identified -- different parameter values give exactly the same distribution.
#>   A support of at least 3 points is required.
```
