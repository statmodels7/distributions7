# Lower-Order Parameter-Scale Derivatives for the Chain Rule

Collects the parameter-scale derivatives of every order strictly below
`order`, in the layout
[`to_link_scale`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
expects.

## Usage

``` r
link_scale_lower_orders(distrib, y, theta, expected, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- expected:

  Logical; whether the derivatives wanted are expected ones.

- order:

  The order being assembled, 1 to 4.

## Value

A list of length `order - 1`, its \\m\\-th element the named list of
order-\\m\\ parameter-scale derivatives.

## Details

Faa di Bruno mixes all lower orders into each link-scale component, so
they must all be to hand before the assembly starts. For **expected**
derivatives the first-order slot is filled with zeros rather than with
the score, because \\\mathbb{E}\[\ell_i\] = 0\\: the term is genuinely
absent, not merely unavailable, which is also why the expected
information transforms as a plain congruence.

## See also

[`to_link_scale`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
