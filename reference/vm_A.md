# The Mean Resultant Length of a von Mises

\\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\, the expected cosine of the
deviation from the location.

## Usage

``` r
vm_A(kappa)
```

## Arguments

- kappa:

  The concentration, a positive numeric vector.

## Value

A numeric vector in \\(0, 1)\\.

## Details

Both Bessel functions are taken exponentially scaled, so the factor
\\e^{\kappa}\\ they share cancels in the ratio and the result stays
finite for a concentration of any size, where the unscaled functions
overflow past about \\\kappa = 700\\.

## See also

[`vonmises_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises_distrib.md),
[`vm_dA`](https://statmodels7.github.io/distributions7/reference/vm_dA.md)
