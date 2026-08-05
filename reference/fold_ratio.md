# The Block Ratios of a Folded Density

Returns a memoized function giving \\d^B L / L\\ for any block, which is
what
[`log_deriv`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
consumes.

## Usage

``` r
fold_ratio(parent, x, theta, order, params, w)
```

## Arguments

- parent:

  The wrapped distribution.

- x:

  A numeric vector of folded observations.

- theta:

  A named list of the parent's parameters.

- order:

  The highest order needed, 1 to 4.

- params:

  The parent's parameter names, in declaration order.

- w:

  The weight of the positive preimage, from
  [`fold_parts`](https://statmodels7.github.io/distributions7/reference/fold_parts.md).

## Value

A function of one block, returning that ratio's vector.

## Details

The ratio is the parent's complete Bell polynomial at each preimage,
weighted by which preimage the point came from: \\w\\(d^B f(x)/f(x)) +
(1-w)\\(d^B f(-x)/f(-x))\\. Both parent evaluations are fetched once and
the result memoized, since a partition sum at fourth order asks for the
same blocks repeatedly.

## See also

[`folded`](https://statmodels7.github.io/distributions7/reference/folded.md),
[`fold_parts`](https://statmodels7.github.io/distributions7/reference/fold_parts.md)
