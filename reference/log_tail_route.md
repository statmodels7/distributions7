# The Log-Tail Route of the CDF Derivatives

Called from the bodies of
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md),
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md),
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
and
[`distrib_deriv4_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
before dispatch. On the log scale it returns the derivatives of \\\log
P\\, \\P\\ the requested tail, from
[`ls_log_tail_derivs()`](https://statmodels7.github.io/distributions7/reference/ls_log_tail_derivs.md)
for the eight location-scale families listed in `log_tail_ls`, and from
[`quad_log_tail_derivs()`](https://statmodels7.github.io/distributions7/reference/quad_log_tail_derivs.md)
for the families listed in `log_tail_quad` (and, for the lower tail,
`log_tail_quad_lower`) at the quantiles where \\\log P\\ is below
`log_tail_cut`, the other quantiles going to the family's own method.
Elsewhere it returns `NULL` and the generic dispatches as usual.

`log_tail_ls` maps each location-scale class to its base and to the map
\\w\\ of its second parameter; `log_tail_quad` lists the classes of the
quadrature route for both tails and `log_tail_quad_lower` those for the
lower tail, which adds the two Weibull families; `log_tail_cut` is the
threshold, \\-5\\; `log_tail_off` holds the switch that turns the route
off while quantiles are sent back to the family's method.

The quantiles are selected on
[`tail_cdf()`](https://statmodels7.github.io/distributions7/reference/tail_probability.md),
the value the family's method computes first. While the method runs on
the quantiles that are not far, `log_tail_stash` holds that value with
the class, quantiles, parameters and tail it belongs to, and
[`tail_cdf()`](https://statmodels7.github.io/distributions7/reference/tail_probability.md)
returns it for a call with the same arguments instead of evaluating the
distribution function again. A parameter of length one is passed to the
method as a scalar.

## Usage

``` r
log_tail_route(distrib, q, theta, order, lower.tail, log, generic)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- q:

  A numeric vector of quantiles, already recycled by the generic.

- theta:

  A named list of parameters, already aligned by the generic.

- order:

  The derivative order, 1 to 4.

- lower.tail:

  A single logical; `TRUE` selects the lower tail and `FALSE` the upper
  tail.

- log:

  A single logical; `TRUE` selects the derivatives of the log
  probability and `FALSE` those of the probability.

- generic:

  The generic of that order, called back on the quantiles that stay with
  the family's method: those that are not far, and those whose tail
  reaches outside the double range, which
  [`quad_log_tail_derivs()`](https://statmodels7.github.io/distributions7/reference/quad_log_tail_derivs.md)
  marks.

## Value

A named list of numeric vectors, keyed as
[`deriv_names(distrib@params, order)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
(in the order of
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order two), or `NULL`.

## Examples

``` r
d <- gaussian1_distrib()
# log S = log(1 - Phi(z)) at z = 40, where 1 - Phi(z) underflows
distributions7:::log_tail_route(d, 40, list(mu = 0, sigma = 1), 1L,
                                FALSE, TRUE, distrib_grad_cdf)
#> $mu
#> [1] 40.02497
#> 
#> $sigma
#> [1] 1600.999
#> 
```
