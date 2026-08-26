# Derivative Components of the Beta-Binomial in Mean and Dispersion

Returns the components of a derivative of order one to four of the
beta-binomial log-mass with respect to \\\mu\\ and \\\sigma\\, by
carrying the shape parametrization's own derivatives through the map
\\(\alpha, \beta) = (\mu/\sigma, (1-\mu)/\sigma)\\.

## Usage

``` r
betabinom1_components(distrib, y, theta, order)
```

## Arguments

- distrib:

  A
  [`BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/BetaBinom1Distrib.md)
  object, read for its `size` and `params`.

- y:

  A numeric vector of counts in \\\\0, \dots, n\\\\.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. `mu` must lie in \\(0, 1)\\ and
  `sigma` be strictly positive.

- order:

  The derivative order, an integer from 1 to 4.

## Value

A named list of component vectors, one per distinct multi-index of the
given order and keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them, so four components at order 3 and five at order 4. Each has
the recycled length of the inputs.

## Details

The shape parametrization
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
carries closed derivatives at every order, each a difference of
polygammas, and this parametrization is that one composed with a map.
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
runs the Faa di Bruno partition sum over the map, with the map's own
partials supplied by
[`md_betabinom1()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
as a keyed table; a key absent from the table is an exact zero. Both
shapes are linear in \\\mu\\ at fixed \\\sigma\\, so every partial
carrying two or more \\\mu\\ vanishes and the sum is short.

The construction is the one
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
runs, used here by a family whose first two orders are written out by
hand in a compiled kernel. A fresh `BetaBinom2Distrib` object is built
on each call, which costs one S7 construction per evaluation.

## See also

[`distrib_deriv3.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BetaBinom1Distrib.md)
and
[`distrib_deriv4.BetaBinom1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BetaBinom1Distrib.md),
which call this;
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
for the partition sum;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
for the family.
