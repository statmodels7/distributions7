# Beta-Binomial Third and Fourth Derivatives in Mean and Dispersion

Closed form at both orders. The shape parametrization carries closed
derivatives at every order, and this one is that one at \\a =
\mu/\sigma\\ and \\b = (1-\mu)/\sigma\\, so the partition sum of
[`chain_derivatives`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
over the map delivers them. Every partial of the map with two or more
\\\mu\\ vanishes, both shapes being linear in \\\mu\\ at fixed
\\\sigma\\.

## Arguments

- distrib:

  A `BetaBinom1Distrib` object.

- y:

  A numeric vector of counts.

- theta:

  A list containing `mu` and `sigma`.

- expected:

  Logical; if `TRUE`, the expected derivatives.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  The approximation used when `expected` is `TRUE`.

- nsim:

  Monte Carlo draws when `approx = "mc"`.

- ...:

  Unused.

## Value

A named list of third-derivative components.

A named list of fourth-derivative components.

## See also

[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
