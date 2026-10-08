# The Scalar Route of a Distribution

Returns how the package's compiled scalar entry points address a
distribution: the name that the C function `d7_scalar_id()` recognizes,
and the constants of the distribution that follow its parameters in the
parameter vector those entry points read. A consumer that resolves the
entry points with `R_GetCCallable()` uses the name to obtain the
distribution's identifier once, and appends the constants to the
parameters of every observation.

## Usage

``` r
distrib_scalar_route(distrib, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from `distrib`.

- ...:

  Passed to methods.

## Value

`NULL` when the registry does not cover the distribution; otherwise a
list with `name`, a single character string, and `constants`, a named
list of numeric vectors, each of length 1 or of the number of
observations, in the order they follow the parameters.

## Details

The entry points are `d7_scalar_id(name)`,
`d7_score_curv(id, k, y, th, out)`, which writes the score and the \\(k,
k)\\ second derivative of the log-density in parameter \\k\\, and
`d7_info_dinfo(id, k, y, th, out)`, which writes the \\(k, k)\\ expected
second derivative and its derivative in the same parameter, all on the
parameter scale, with `k` counted from zero over `distrib@params` and
`th` the parameters of one observation followed by its constants.
`d7_logpdf(id, y, th)` returns the log-density of one observation.
`d7_scalar_thread_safe(id)` returns 1 when the entries of a distribution
never reach the R API, so that they may be called from a worker thread,
and 0 otherwise.

The default method returns the class name when the registry covers the
class and no constants; the binomial and the two beta-binomials carry
their size as a constant.

A wrapped family is named `"<wrapper class>[:aux]|<inner name>"`, `aux`
being what the wrapper needs besides its constants and the inner name
being the route's name of the wrapped family, which may itself be a
wrapper; its constants are the inner family's followed by the wrapper's.
For
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
`aux` is the mask of the fixed parameters (bit `j - 1` for the `j`-th
parameter of the inner family) and the wrapper's constants are the fixed
values in the inner family's order. For
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
there is no `aux` and no constant of the wrapper's own, the probability
being the last parameter. For
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md),
`aux` is the code of the transformer, which must be one of the twelve
ready-made ones, and the wrapper's constants are the transformer's
parameters; a transformer built with
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
has no route.

## Examples

``` r
distrib_scalar_route(gaussian1_distrib())
#> $name
#> [1] "Gaussian1Distrib"
#> 
#> $constants
#> list()
#> 
distrib_scalar_route(binomial_distrib(size = 10))
#> $name
#> [1] "BinomialDistrib"
#> 
#> $constants
#> $constants$size
#> [1] 10
#> 
#> 
is.null(distrib_scalar_route(mvgaussian1_distrib(n_dim = 2)))
#> [1] TRUE
```
