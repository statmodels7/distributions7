# The Derivative of the Mean Resultant Length

\\A'(\kappa) = 1 - A(\kappa)/\kappa - A(\kappa)^2\\, which is also the
variance of \\\cos(Y - \mu)\\ and therefore positive.

## Usage

``` r
vm_dA(kappa, A = vm_A(kappa))
```

## Arguments

- kappa:

  The concentration, a positive numeric vector.

- A:

  The value of
  [`vm_A`](https://statmodels7.github.io/distributions7/reference/vm_A.md)
  at `kappa`, passed in when it has already been computed.

## Value

A numeric vector.

## Details

The identity follows from \\I_0' = I_1\\ and \\I_1' = I_0 -
I_1/\kappa\\, so no further Bessel evaluation is needed.

## See also

[`vonmises_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises_distrib.md),
[`vm_A`](https://statmodels7.github.io/distributions7/reference/vm_A.md)
