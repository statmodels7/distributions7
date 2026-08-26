# Second-Response Mixed Derivatives Through a Map

Carries the parent's
[`distrib_cross2_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.md)
block onto a new parametrization by the first-order chain rule, \\\sum_i
(\partial\theta_i/\partial\alpha_a)\\
\partial\ell^{(yy)}/\partial\theta_i\\. It is the same expansion
[`mapped_cross_y()`](https://statmodels7.github.io/distributions7/reference/mapped_cross_y.md)
takes, read one derivative further in the response.

## Usage

``` r
mapped_cross2_y(distrib, parent, th_par, maps, y)
```

## Arguments

- distrib:

  The distribution in the new parametrization. Only its `params` are
  read, to key the result.

- parent:

  The parent distribution, which supplies the block.

- th_par:

  A named list of the parent's parameters, evaluated at the new ones.

- maps:

  The map's keyed partial tables, from
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md).

- y:

  A numeric vector of observations.

## Value

A named list with one numeric vector per new parameter, keyed by
`distrib@params`.

## Details

Only the map's FIRST partials enter, so the second-partial tables
[`mapped_theta2()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md)
needs are never touched. A missing key is an exact zero and is skipped.

## Notation

\\\ell\\ is the log-density of one observation, \\\theta_i\\ a parameter
of the PARENT and \\\alpha_a\\ one of the new parametrization, so that
\\\theta = \theta(\alpha)\\ is the map. \\\ell^{(y)}\\ and
\\\ell^{(yy)}\\ are the first and second derivatives in the response.

## See also

[`mapped_theta2()`](https://statmodels7.github.io/distributions7/reference/mapped_theta2.md)
for the second order in \\\theta\\,
[`distrib_cross2_y.ReparamContinuousDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cross2_y.ReparamContinuousDistrib.md),
its caller, and
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
for the partial tables.

## Examples

``` r
# gaussian2 carries (mu, sigma2) against a parent in (mu, sigma), so the
# scale component picks up d(sigma)/d(sigma2) = 1 / (2 sqrt(sigma2)).
d <- gaussian2_distrib()
y <- c(-0.7, 0.3, 1.4)
theta <- list(mu = 0.3, sigma2 = 1.44)

vapply(distrib_cross2_y(d, y, theta), function(z) z[1], numeric(1))
#>        mu    sigma2 
#> 0.0000000 0.4822531 

# The parent's own block, carried by that one factor.
par_block <- distrib_cross2_y(gaussian1_distrib(), y,
                              list(mu = 0.3, sigma = 1.2))
c(mu = par_block$mu[1], sigma2 = par_block$sigma[1] / (2 * sqrt(1.44)))
#>        mu    sigma2 
#> 0.0000000 0.4822531 
```
