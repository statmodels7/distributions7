# Zero-Adjusted Distribution Object

Creates a zero-adjusted version of an existing distribution: the
probability of a zero becomes a parameter of its own, \\\pi\\ (parameter
`za`), and everything else is left to the parent. What that means
depends on the parent's type, and the constructor dispatches on it:

- **Discrete** (support including 0): a **hurdle** model. The mass the
  parent puts at zero is removed, the parent is renormalised over the
  positive values, and \\\pi\\ takes its place.

- **Continuous**: a **mixed** distribution. Nothing has to be removed,
  since \\P(Y = 0) = 0\\ already; a point mass \\\pi\\ is placed at zero
  and the density is scaled by \\1-\pi\\.

Zero-adjustment is the right wrapper when zeros come from their own
mechanism — no claim filed, no purchase made, no rainfall — and the
parent describes only what happens once that mechanism has been passed.
When the zeros instead come partly from the parent itself and partly
from a separate one, so that no single zero can be attributed, use
[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md).

## Usage

``` r
zero_adjusted(distrib, link_za = logit_link())
```

## Arguments

- distrib:

  An object inheriting from `discrete_distrib` (with 0 in its support)
  or from `continuous_distrib`.

- link_za:

  A link function object for the zero probability \\\pi\\. Defaults to
  [`logit_link`](https://rdrr.io/pkg/linkfunctions7/man/logit_link.html).

## Value

An S7 object of class `ZeroAdjustedDiscreteDistrib` or
`ZeroAdjustedContinuousDistrib`.

## Details

**Discrete parent (hurdle).** \$\$ P(Y=y; \theta, \pi) = \begin{cases}
\pi & y = 0 \\ (1 - \pi)\dfrac{f(y; \theta)}{1 - f(0; \theta)} & y \> 0
\end{cases} \$\$ The division by \\1 - f(0;\theta)\\ is the truncation:
it is what distinguishes this from zero-inflation, and what makes
\\\theta\\ the parameters of a law on the positive integers rather than
of the original count process. The log-likelihood separates into a
Bernoulli part in \\\pi\\ and a truncated part in \\\theta\\, so the
mixed blocks of the information matrix are exactly zero and the two
halves could in principle be fitted separately.

**Continuous parent (mixed).** \$\$f_Y(0) = \pi, \qquad f_Y(y) = (1-\pi)
f_W(y;\theta) \quad (y \neq 0)\$\$ Here \\f_Y\\ is a density with
respect to Lebesgue measure plus a point mass at zero, so it integrates
to \\1 - \pi\\: the remainder is the atom, which
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
reports and
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
accounts for. The classical members of this family — zero-adjusted
gamma, inverse Gaussian or lognormal for semicontinuous data,
zero-adjusted beta for proportions — all have a parent supported on
\\(0, \infty)\\ or \\(0,1)\\, where zero sits on the boundary. A parent
supported on the whole line is accepted too and gives a "spike at zero"
model, but note that \\y = 0\\ then no longer identifies its own
mechanism: the atom sits where the density is positive.

**Choosing between the two wrappers.** Zero-inflation can only add
zeros, since \\P(Y=0) = \zeta + (1-\zeta)f(0) \> f(0)\\; the hurdle
replaces \\f(0)\\ outright and so also covers the case of *fewer* zeros
than the parent implies. Where both apply they are not nested, and they
differ in interpretation more than in fit: zero-inflation keeps
\\\theta\\ as the parameters of the original count process, while the
hurdle re-reads them as those of a truncated one. Prefer the hurdle when
a zero is observable evidence of a distinct decision, and zero-inflation
when the two kinds of zero are genuinely indistinguishable.

**What the parent must be.** A discrete parent must have 0 in its
support: with \\f(0) = 0\\ there is no mass to remove. Construction also
fails when the result would not be identified:

- the parent already models a probability of zero — zero-truncating a
  zero-inflated or zero-adjusted distribution cancels its zero parameter
  out of the likelihood entirely, leaving an identically zero score;

- the support is too small to carry one more parameter: a distribution
  on \\k\\ points has \\k-1\\ free probabilities, so at least
  `n_params + 2` support points are needed. Zero-adjusting a Bernoulli
  leaves the truncated part concentrated on \\\\1\\\\, and `mu` vanishes
  from the likelihood.

A continuous parent whose support does not reach zero (say \\(2, 5)\\)
is accepted with a warning: the atom is then disconnected from the rest
of the distribution, which is legitimate but rarely intended.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
for the mixture counterpart,
[`distrib_atoms`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md),
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md).

## Examples

``` r
# Hurdle Poisson: the mass at zero is exactly za, not dpois(0, mu)
zap <- zero_adjusted(poisson_distrib())
distrib_pdf(zap, 0:5, list(mu = 3, za = 0.3))
#> [1] 0.3000000 0.1100310 0.1650464 0.1650464 0.1237848 0.0742709

# Semicontinuous data: a spike at zero and a gamma above it
zagamma <- zero_adjusted(gamma_distrib())
distrib_atoms(zagamma, list(mu = 2, sigma2 = 1, za = 0.3))
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3
#> 

# The truncated part of a zero-adjusted Bernoulli has no free parameter
try(zero_adjusted(bernoulli_distrib()))
#> Error : zero_adjusted() cannot wrap 'bernoulli': its support has 2 points, so the family has 1 free
#>   probability, while the wrapped distribution would have 2 parameters. They are not
#>   identified -- different parameter values give exactly the same distribution.
#>   A support of at least 3 points is required.
```
