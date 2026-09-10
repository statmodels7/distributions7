# Reject Centered Parameters the Map Cannot Carry

Checks that the direct parameters
[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
has produced are usable, and raises in the centered family's own terms
when they are not.

## Usage

``` r
sn2_reject_unmappable(dp, theta)
```

## Arguments

- dp:

  The direct parameters, as returned by
  [`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md).

- theta:

  The centered parameters the caller supplied, ordered as
  `c("mu", "sigma", "gamma1")`.

## Value

Invisibly `NULL`; raises an error naming `"skew normal2"` if any mapped
parameter is not finite or the mapped scale is not positive.

## Details

Every probability function of
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
evaluates the parent at the mapped parameters, and the parent validates
what it is handed against its own domains. Without this check a caller
who wrote `gamma1` reads an error about `alpha`, a parameter the model
does not have, and about `"skew normal1"`, a family the call does not
name. The message here reports the centered parameter responsible and
the value it holds.

The check states a property of the delegation rather than guarding a
value the public surface can reach.
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
bounds `gamma1` and every generic validates it before dispatch, so a
skewness outside its domain is reported before the map runs; what
remains reachable is a `sigma` and a `gamma1` each inside its own domain
whose implied scale or location leaves the doubles.

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
#>   'gamma1' = 0.9952717 is at the ceiling 0.9952717 of the skewness a skew normal
#>   can carry, where the shape it maps to is not finite.

# A scale and a skewness each inside its own domain whose implied scale
# is not.
th2 <- list(mu = 0, sigma = 1.1e308, gamma1 = 0.99)
try(distributions7:::sn2_theta(th2))
#> Error : Invalid parameter value(s) for the 'skew normal2' distribution:
#>   'sigma' = 1.1e+308 and 'gamma1' = 0.99 imply a scale that is not finite.
```
