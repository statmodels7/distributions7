# Register the Higher CDF Orders on a Mapped Family

Turns a parent and a map into the two or three S7 methods that a mapped
family needs, so that the family states its map once instead of once per
order. Five families are registered through it: the two further Gaussian
parametrizations, the second inverse Gaussian, the second Laplace and
the lognormal.

## Usage

``` r
register_mapped_cdf_k(
  cls,
  parent_fn,
  th_fn,
  md_fn,
  q_fn = identity,
  q_ok = NULL,
  orders = 3:4
)
```

## Arguments

- cls:

  The S7 class to register on.

- parent_fn:

  A function of no arguments returning the parent distribution. The
  parent is given as a function rather than as the object, so that it is
  built at call time; at load time the class that it names may not exist
  yet.

- th_fn:

  A function of `theta` returning the parent's parameters.

- md_fn:

  A function of `theta` and the order, returning the map's keyed partial
  tables up to that order.

- q_fn:

  The transformation of the response, for a parent that is the same law
  on a transformed scale. The default is `identity`, and the lognormal
  uses `log`.

- q_ok:

  `NULL`, the default, or a function of `q` returning a logical vector
  that is `FALSE` where `q_fn` is undefined and the distribution
  function is identically zero, as at \\q \le 0\\ for the lognormal.
  There every derivative is set to zero, and `q_fn` is evaluated at 1.

- orders:

  An integer vector of the orders to register, `3:4` by default. A
  family whose written-out route stops below the fourth order takes the
  rest here.

## Value

Invisibly `NULL`. The function is called for the registration.

## Details

A family registers more than the top two orders by setting `orders`. The
second inverse Gaussian is registered at 2 to 4 because its written-out
route in `cdf_derivatives_families.R` registers only the gradient, so
its Hessian is taken from the chain here.

`force(o)` inside the factory prevents the registrations from sharing
one order.

## See also

[`mapped_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/mapped_cdf_deriv_k.md),
the body that it registers;
[`distrib_deriv3_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3_cdf.md)
for the generics.
