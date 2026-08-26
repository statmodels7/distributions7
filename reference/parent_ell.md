# Look Up the Parent's Derivative Components by Block

Fetches the parent's derivatives up to `max_order` and returns a
function giving the component belonging to any block.

## Usage

``` r
parent_ell(parent, y, theta, max_order, params)
```

## Arguments

- parent:

  The parent distribution.

- y:

  A numeric vector of observations.

- theta:

  A named list of the parent's parameters.

- max_order:

  The highest order needed, 1 to 4.

- params:

  The parent's parameter names, in declaration order.

## Value

A function of one block, returning that component's vector.

## Details

The orders are fetched once per call rather than once per block. A
partition sum at fourth order asks for the same handful of components
repeatedly, and the parent's derivative may itself be expensive.

## See also

[`canon_key()`](https://statmodels7.github.io/distributions7/reference/canon_key.md)
