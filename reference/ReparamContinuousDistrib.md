# S7 Classes for a Reparametrized Distribution

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
returns objects of these classes. Each of them represents the same law
as its parent, written in different coordinates. There is one class per
kind of parent, so that a continuous parent keeps the defaults
registered on
[`continuous_distrib()`](https://statmodels7.github.io/distributions7/reference/continuous_distrib.md)
and a discrete one those of
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).
The two classes differ only in their parent, and the methods of this
file are registered on both from the same body, except
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
which only the continuous class registers.

## Usage

``` r
ReparamContinuousDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  parent_distrib = NULL,
  reparam_map = function() NULL,
  reparam_derivs = function() NULL
)

ReparamDiscreteDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  parent_distrib = NULL,
  reparam_map = function() NULL,
  reparam_derivs = function() NULL
)
```

## Arguments

- distrib_name:

  The name of the family, a single string.

- dimension:

  `"univariate"` or `"multivariate"`. It is always `"univariate"` for a
  reparametrization, since
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
  signals an error for a multivariate parent.

- bounds:

  A length-two numeric vector, the support.

- params:

  A character vector naming the new parameters.

- params_interpretation:

  A named character vector describing each new parameter.

- n_params:

  The number of new parameters.

- params_bounds:

  A named list of length-two numeric vectors, each of which gives the
  open interval in which the corresponding new parameter lies.

- link_params:

  A named list of linkfunctions7 links, one per new parameter. The name
  is `link_params` on the object and `links` in
  [`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)'s
  signature.

- params_smooth:

  A named logical vector indicating the parameters with respect to which
  the log-density is differentiable.

- parent_distrib:

  The distribution that is rewritten.

- reparam_map:

  The map from the new parameters to the parent's, a function of one
  named list returning another.

- reparam_derivs:

  The function returning the map's keyed partial tables, as
  [`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md)
  consumes them.

## Value

An object of class `ReparamContinuousDistrib` or
`ReparamDiscreteDistrib`, carrying every property of a `distrib` plus
the three above.

## Details

The methods read three properties beyond a distribution's own.
`parent_distrib` is the law that is rewritten, and the methods delegate
to it for the density, the distribution function, the quantile function,
the generator and the response derivatives, all of which a change of
parametrization leaves unchanged. `reparam_map` carries the new
parameters to the parent's, and `reparam_derivs` carries the map's
partial derivatives, on which the exactness of the parameter derivatives
depends.

## See also

[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
the constructor that fills these in;
[`is_reparam()`](https://statmodels7.github.io/distributions7/reference/is_reparam.md)
to test for either class;
[`reparam_theta()`](https://statmodels7.github.io/distributions7/reference/reparam_theta.md)
and
[`reparam_tables()`](https://statmodels7.github.io/distributions7/reference/reparam_tables.md),
the two readers.

## Examples

``` r
# The class below is the one that a continuous parent gives.
d <- reparametrize(
  gaussian1_distrib(),
  map = function(psi) list(mu = psi$mu, sigma = sqrt(psi$sigma2)),
  params = c("mu", "sigma2"),
  bounds = list(mu = c(-Inf, Inf), sigma2 = c(0, Inf)),
  links = list(mu = linkfunctions7::identity_link(),
               sigma2 = linkfunctions7::log_link())
)
class(d)
#> [1] "distributions7::ReparamContinuousDistrib"
#> [2] "distributions7::continuous_distrib"      
#> [3] "distributions7::distrib"                 
#> [4] "S7_object"                               
d@parent_distrib@distrib_name
#> [1] "gaussian1"
```
