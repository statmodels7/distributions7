# The Order of Differentiability of a Family in Each Parameter

One integer per parameter: the largest `k` such that the log-density has
`k` continuous derivatives in that parameter everywhere in `y`. `Inf`
says the parameter is smooth, which is the answer for every family that
declares no
[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md).

## Usage

``` r
params_order(distrib, ...)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- ...:

  Passed to
  [`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md).

## Value

A named numeric vector, one entry per element of `distrib@params`.

## Details

The order is DEDUCED and not declared. Where a family declares a
decomposition \\c(\theta)\phi(v(y,\theta))\\, a parameter the kink moves
with – one whose \\\partial v/\partial \theta_p\\ is not identically
zero – inherits the order of \\\phi\\, and every other parameter stays
`Inf`. Whether \\\partial v/\partial \theta_p\\ vanishes is measured at
probe values rather than assumed, so a Huber cut-off held at zero
reports the scale as smooth, which it then is.

A family that declares no decomposition but records a parameter as
non-smooth through `params_smooth` reports `NA` for it: there is a kink
and its order has not been established. A consumer must treat that as
unusable rather than as smooth, which is what
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
and the fitting layer do. This is the state of every wrapper of a kinked
family today.

The consumers and the order each one needs:

|  |  |
|----|----|
| 1 | a score continuous in the coefficients |
| 2 | IWLS, Newton, [`vcov()`](https://rdrr.io/r/stats/vcov.html), and the determinant of a Laplace criterion |
| 3 | the exact gradient of a marginal criterion |
| 4 | its exact Hessian |
| 5 | the fourth derivative of a filtered predictor |

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md),
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md),
[`param_smoothness()`](https://statmodels7.github.io/distributions7/reference/param_smoothness.md).

## Examples

``` r
# smooth everywhere
params_order(gaussian1_distrib())
#>    mu sigma 
#>   Inf   Inf 

# the Laplace: the location carries the kink, the scale does not
params_order(laplace_distrib())
#>    mu sigma 
#>     0   Inf 
```
