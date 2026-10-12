# Split a Wrapper's Parameter From Its Parent's

Separates a full `theta` into the parent distribution's parameters and
the single mixing parameter the wrapper adds. Both zero wrappers append
their parameter last, so the split is positional and uses no parameter
names. The constructors reject a parent that already has a parameter
with the wrapper's name (`zi` or `za`).

## Usage

``` r
split_mix_theta(distrib, theta)
```

## Arguments

- distrib:

  A `ZeroInflatedDistrib` or a zero-adjusted object, whose `n_params`
  gives the position of the wrapper's parameter.

- theta:

  A named list of parameters, already aligned, the parent's followed by
  the wrapper's.

## Value

A named list with `orig`, the parent's parameters as a named list, and
`mix`, the wrapper's parameter as a numeric vector.

## See also

[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md),
whose methods all begin with this.

## Examples

``` r
d <- zero_inflated(poisson_distrib())
theta <- list(mu = 3, zi = 0.25)
p <- distributions7:::split_mix_theta(
  d, distributions7:::align_theta(d, theta))
str(p)
#> List of 2
#>  $ orig:List of 1
#>   ..$ mu: num 3
#>  $ mix : num 0.25

# The first block is exactly what the parent expects.
all.equal(distrib_pdf(d@parent_distrib, 0:3, p$orig), dpois(0:3, 3))
#> [1] TRUE
```
