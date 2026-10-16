# Reject Centered Parameters That the Map Cannot Carry

Checks that the direct parameters produced by
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
are usable, and signals an error in the centered family's own terms when
they are not.

## Usage

``` r
sn2_reject_unmappable(dp, theta)
```

## Arguments

- dp:

  The direct parameters, as returned by
  [`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md).

- theta:

  The centered parameters that the caller supplied, ordered as
  `c("mu", "sigma", "gamma1")`.

## Value

Invisibly `NULL`; raises an error naming `"skew normal2"` if any mapped
parameter is not finite or the mapped scale is not positive.

## Details

Every probability function of
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
evaluates the parent at the mapped parameters, and the parent validates
what it is handed against its own domains. Without this check a caller
who wrote `gamma1` reads an error about `alpha`, a parameter that the
model does not have, and about `"skew normal1"`, a family that the call
does not name. The message reports the centered parameter that is
responsible and the value that it holds.

The check covers the failures of the delegation itself.
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
bounds `gamma1`, and every generic checks the parameters against their
domains before it calls the method, so a skewness outside its domain is
reported before the map runs. The remaining case is a `sigma` and a
`gamma1`, each inside its own domain, whose implied scale or location is
not finite.

## See also

[`sn2_theta()`](https://statmodels7.github.io/distributions7/reference/sn2_theta.md),
which calls it, and
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
for the map.

## Examples

``` r
# A skewness at the ceiling is reported in the centered family's own terms.
th <- list(mu = 0, sigma = 1, gamma1 = distributions7:::sn_max_skew())
try(distributions7:::sn2_theta(th))
#> Error : Invalid parameter value(s) for the 'skew normal2' distribution:
#>   'gamma1' = 0.9952717 is at the ceiling 0.9952717 of the skewness that a skew normal
#>   can carry, where the shape to which it maps is not finite.

# A scale and a skewness that each lie inside their own domains can still
# imply a scale that is not finite.
th2 <- list(mu = 0, sigma = 1.1e308, gamma1 = 0.99)
try(distributions7:::sn2_theta(th2))
#> Error : Invalid parameter value(s) for the 'skew normal2' distribution:
#>   'sigma' = 1.1e+308 and 'gamma1' = 0.99 imply a scale that is not finite.
```
