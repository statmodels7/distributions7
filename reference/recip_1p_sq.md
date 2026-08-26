# Derivatives of the Reciprocal of One Plus a Square

Returns \\\partial^j (1+a^2)^{-1}/\partial a^j\\ for \\j = 0, \ldots,
3\\: \\(1+a^2)^{-1}\\, \\-2a(1+a^2)^{-2}\\, \\(6a^2-2)(1+a^2)^{-3}\\ and
\\24a(1-a^2)(1+a^2)^{-4}\\. This is the factor the skew normal's shape
derivative carries, and writing it out saves a Faa di Bruno pass over a
quantity whose derivatives are four short expressions.

## Usage

``` r
recip_1p_sq(a, j)
```

## Arguments

- a:

  A numeric vector, the shape.

- j:

  The order, 0 to 3. Any other value returns `NULL`,
  [`switch()`](https://rdrr.io/r/base/switch.html) falling through;
  three is the highest the fourth-order cdf derivative needs.

## Value

A numeric vector the length of `a`.

## See also

[`sn_cdf_std_derivs()`](https://statmodels7.github.io/distributions7/reference/sn_cdf_std_derivs.md),
the one consumer.

## Examples

``` r
# The first derivative at a = 2 is -2a / (1 + a^2)^2 = -0.16.
distributions7:::recip_1p_sq(2, 1)
#> [1] -0.16
```
