# Lower-Order Parameter-Scale Derivatives for the Chain Rule

Collects the parameter-scale derivatives of every order strictly below
`order`, in the layout that
[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
expects.

## Usage

``` r
link_scale_lower_orders(distrib, y, theta, expected, order, method = NULL)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- expected:

  Logical; whether the derivatives to compute are expected ones.

- order:

  The order to assemble, 1 to 6. Order 6 is observed only, and
  `expected = TRUE` with it signals an error.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
  passed to the expected derivatives of the lower orders so that every
  order of a call is obtained by the same method. It is read only when
  `expected` is `TRUE`.

## Value

A list of length `order - 1`, in which the \\m\\-th element is the named
list of order-\\m\\ parameter-scale derivatives.

## Details

Faa di Bruno mixes all lower orders into each link-scale component, so
they must all be to hand before the assembly starts. For expected
derivatives the first-order slot is filled with zeros rather than with
the score, because \\\mathbb{E}\[\ell_i\] = 0\\. The term is absent from
the expected derivatives, which is also why the expected information
transforms as a plain congruence.

## See also

[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
