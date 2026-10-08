# Generalized Gamma Derivatives in the Response

Return the derivatives of the log-density in the response to order four,
and the mixed derivatives of orders one and two in the response and one
and two in \\(a, d, p)\\, each a closed form from its own compiled
kernel, written as
[`distrib_deriv3.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
describes.

## Arguments

- distrib:

  A `GenGamma1Distrib` object, from
  [`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `y`, all strictly positive.

- scale:

  `"parameter"` or `"link"`, for the mixed derivatives; the link scale
  is applied by the generic.

- ...:

  Unused.

## Value

A numeric vector for the derivatives in the response; a named list,
keyed by parameter or by parameter pair in
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)'s
order, for the mixed derivatives.

## See also

[`distrib_gradient.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma1Distrib.md)
for the derivatives in the parameters.

## Examples

``` r
d <- gengamma1_distrib()
th <- list(a = 2, d = 3, p = 1.5)
y <- c(0.5, 1.5, 4)

# The score in y written out: ((d - 1) - p (y/a)^p) / y.
all.equal(distrib_grad_y(d, y, th), (2 - 1.5 * (y / 2)^1.5) / y)
#> [1] TRUE

distrib_hess_y_hess(d, y, th)$p_p
#> [1]  1.0519092 -0.2630819 -0.9073821
```
