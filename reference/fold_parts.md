# The Two Preimages of a Folded Point

The parent's density at \\+x\\ and at \\-x\\, their sum, and the weight
\\w = f(x)/L\\ the first carries.

## Usage

``` r
fold_parts(parent, x, theta)
```

## Arguments

- parent:

  The wrapped distribution.

- x:

  A numeric vector of folded observations.

- theta:

  A named list of the parent's parameters.

## Value

A list with `fp`, `fm`, `L` and `w`.

## Details

Every method of
[`FoldedDistrib`](https://statmodels7.github.io/distributions7/reference/FoldedDistrib.md)
needs the same four quantities, and computing them once keeps the
parent's density from being evaluated twice per call. Points outside the
folded support contribute nothing and are returned with a zero density
rather than being dropped, so that the result aligns with the input.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md)
