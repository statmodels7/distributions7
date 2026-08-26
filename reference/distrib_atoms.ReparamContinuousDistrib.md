# Atoms of a Reparametrized Distribution

The parent's, unchanged: a reparametrization does not move mass.
Delegating matters because
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
consults this to decide which of its checks apply, and a reparametrized
mixed family would otherwise be validated as though it were purely
continuous.

## Usage

``` r
reparam_atoms(distrib, theta, ...)
```

## Arguments

- distrib:

  A reparametrized distribution, from
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

- theta:

  A named list of the new parameters, on the new parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list with components `y` (the atom locations) and `p` (their
probabilities), both numeric vectors of the same length, possibly of
length zero.

## See also

[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for the generic and what consults it;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md);
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

## Examples

``` r
# A reparametrized zero-adjusted gamma keeps its atom at zero.
za <- zero_adjusted(gamma2_distrib())
d <- reparametrize(
  za,
  map = function(psi) list(mu = psi$mu, sigma2 = psi$s, za = psi$p0),
  params = c("mu", "s", "p0"),
  bounds = list(mu = c(0, Inf), s = c(0, Inf), p0 = c(0, 1)),
  links = list(mu = linkfunctions7::log_link(),
               s = linkfunctions7::log_link(),
               p0 = linkfunctions7::logit_link())
)
distrib_atoms(d, list(mu = 2, s = 1, p0 = 0.3))
#> $y
#> [1] 0
#> 
#> $p
#> [1] 0.3
#> 
```
