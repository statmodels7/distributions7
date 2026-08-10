# Register the Two New Orders on a Mapped Family

Turns the parent and the map into the two methods, so that a family
states its map once instead of twice.

## Usage

``` r
register_mapped_cdf_k(
  cls,
  parent_fn,
  th_fn,
  md_fn,
  q_fn = identity,
  orders = 3:4
)
```

## Arguments

- cls:

  The S7 class.

- parent_fn:

  A function of no arguments returning the parent.

- th_fn:

  A function of `theta` returning the parent's parameters.

- md_fn:

  The map's table function.

- q_fn:

  The transformation of the response, when the parent is the same law on
  a transformed scale. The identity by default.

- orders:

  The orders to register, 3 and 4 by default. A family whose written-out
  route stops below the fourth order takes the rest here.

## Value

Invisibly `NULL`; called for the registration.
