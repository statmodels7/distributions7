# Random Parameters Inside a Family's Own Domain

Draws one parameter value per entry of `distrib@params`, each inside
that parameter's own interval in `params_bounds` and strictly away from
its ends. A bounded interval is sampled across its width, a half-line is
sampled from a positive draw offset from the finite end, and an
unbounded parameter is sampled around zero.

The draws ignore the response entirely, so they are a fallback and never
an estimate: being of order one whatever the data, they suit a shape or
a probability and are the wrong order for a location or a scale on a
response of any size.
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
uses them for the starting values after the first, and
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
uses them to probe a family at a parameter nobody chose.

## Arguments

- distrib:

  An object inheriting from `distrib`, read for `params` and
  `params_bounds`.

- ...:

  Unused.

## Value

A named list with one numeric value per parameter, named and ordered as
`distrib@params`, every value strictly inside its own bounds.

## See also

[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md),
which draws these for the starts after the first;
[`start_from_moments()`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md)
for the data-based one;
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md),
which probes a family at a random parameter.

## Examples

``` r
set.seed(1)
unlist(generate_random_theta(gaussian1_distrib()))
#>        mu     sigma 
#> -1.406948  2.295531 

# Every draw is admissible by construction: a scale is positive, a
# probability is inside (0, 1).
set.seed(2)
d <- beta1_distrib()
th <- generate_random_theta(d)
rbind(draw = unlist(th), lower = vapply(d@params_bounds, `[`, numeric(1), 1),
      upper = vapply(d@params_bounds, `[`, numeric(1), 2))
#>              mu      phi
#> draw  0.2479058 4.244007
#> lower 0.0000000 0.000000
#> upper 1.0000000      Inf
```
