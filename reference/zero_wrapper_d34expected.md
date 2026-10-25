# The Third and Fourth Derivatives of the Zero Wrappers' Expected Information

Compute the third or fourth derivative of the expected Hessian of
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
families, on the parameter scale, by differentiating the expressions of
[`zero_wrapper_dexpected()`](https://statmodels7.github.io/distributions7/reference/zero_wrapper_dexpected.md)
along each multiset of parameters. The parent enters through its
expected Hessian and its derivatives up to the same order, and, for a
discrete parent, through the log-mass at zero and its derivatives to
order five.

## Usage

``` r
zero_wrapper_d34expected(distrib, y, theta, order, threads = 1L)
```

## Arguments

- distrib:

  A family from
  [`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
  or
  [`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md).

- y:

  The response, read for its length.

- theta:

  The parameters, the parent's followed by the mixing probability.

- order:

  `3L` or `4L`.

- threads:

  The thread count passed to the parent's methods.

## Value

A named list keyed as
[`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
or
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md),
each component of length `length(y)`, on the parameter scale.

## Details

The entries are those of
[`zero_wrapper_dexpected()`](https://statmodels7.github.io/distributions7/reference/zero_wrapper_dexpected.md).
With \\f\\ the parent's mass at zero, \\s_i\\ its score there,
\\E\_{ij}\\ the parent's expected Hessian, \\\zeta\\ the probability and
\\L = \zeta + (1-\zeta) f\\, the zero-inflated entries are
\\(1-\zeta)E\_{ij} + \zeta(1-\zeta) f s_i s_j / L\\, \\-f s_i / L\\ and
\\-(1-f)^2/L - (1-f)/(1-\zeta)\\, and the zero-adjusted discrete ones
are \\(1-\pi)\[E\_{ij}/q + f s_i s_j/q^2\]\\ with \\q = 1-f\\. The
derivative of \\f\\ along a multiset of parameters is formed from the
log-mass derivatives by \\\partial\_{aR} f = \sum\_{U \subseteq R}
\partial_U f\\ \partial\_{a R\setminus U}\log f\\, products by Leibniz's
rule and reciprocals from \\A\\(1/A) = 1\\, through the operations of
[`zw34_memo()`](https://statmodels7.github.io/distributions7/reference/zw34_memo.md).
The orders of the parent's expected Hessian derivatives are read from
its own methods, closed form where the parent has them and the default
stencil otherwise. The order-3 call reads the parent to order 3 and the
order-4 call to order 4.

## See also

[`zero_wrapper_dexpected()`](https://statmodels7.github.io/distributions7/reference/zero_wrapper_dexpected.md),
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md)
