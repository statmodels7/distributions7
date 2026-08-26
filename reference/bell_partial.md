# Partial Bell Polynomials for Orders up to Four

The partial (incomplete) Bell polynomial \\B\_{m,j}\\ evaluated at the
derivatives of the inverse link, for \\m \le 4\\.

## Usage

``` r
bell_partial(m, j, h)
```

## Arguments

- m:

  The total order, 1 to 4.

- j:

  The number of blocks, 1 to `m`.

- h:

  A list with `h[[k]]` the \\k\\-th derivative of the inverse link
  evaluated at \\\eta\\, as a numeric vector.

## Value

A numeric vector, the polynomial evaluated element-wise.

## Details

These are the coefficients Faa di Bruno's formula needs. Because each
parameter carries its own link, the Jacobian of \\\theta \mapsto \eta\\
is diagonal and the multivariate formula factorizes into a product of
univariate ones, so only \\B\_{m,j}\\ for a single variable is required.
They are written out rather than generated: there are ten of them below
order five, and a table cannot be slower or wrong in a way a recursion
could.

## See also

[`link_scale_derivatives()`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md),
[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
