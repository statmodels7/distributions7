# von Mises Distribution Function in the Resultant Length

The concentration parametrization's series, read at the concentration
this one's resultant length implies.

## Arguments

- distrib:

  A `VonMises2Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A named list with `mu` and `rho`.

- ...:

  Unused.

## Value

The distribution function at `q`.

## Details

The map touches the second parameter only and the response not at all,
so the distribution function is the other family's at \\\kappa =
A^{-1}(\rho)\\. What it replaces is the base class's quadrature, one per
observation; see
[`vm_cdf`](https://statmodels7.github.io/distributions7/reference/vm_cdf.md)
for the series and for how many terms it takes.

## See also

[`vm_cdf`](https://statmodels7.github.io/distributions7/reference/vm_cdf.md),
[`vonmises2_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
