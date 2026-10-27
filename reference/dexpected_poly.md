# Polynomial Algebra of Observed Derivatives

Helpers for the derivatives of the expected-information matrix under an
explicit method. A polynomial is a list of terms, each a list with a
numeric `coef` and `f`, a list of sorted integer index tuples, one per
factor \\\ell_T\\.

## Usage

``` r
poly_base(identity, a, b)

poly_d(poly, cc, moved)

poly_order(poly)

poly_eval(poly, getter)
```

## Arguments

- identity:

  `"hessian"` or `"score"`.

- a, b:

  Parameter indices of the pair.

- poly:

  A polynomial.

- cc:

  The index of the differentiating parameter.

- moved:

  Whether the density moves with the parameter.

- getter:

  A function of a sorted index tuple returning the observed derivative,
  a numeric vector.

## Value

`poly_base()` and `poly_d()` return a polynomial, `poly_order()` an
integer and `poly_eval()` a numeric vector.

## Details

`poly_base()` returns the integrand of the identity for the pair
\\(a,b)\\: \\\ell\_{ab}\\ for `"hessian"` and \\-\ell_a\ell_b\\ for
`"score"`. `poly_d()` applies the derivative with respect to the
parameter `cc`: the Leibniz derivative, which replaces each factor
\\\ell_T\\ by \\\ell\_{T\cup c}\\, and, when `moved` is `TRUE`, also the
term \\f\\\ell_c\\ that comes from the density. `poly_order()` is the
largest derivative order that appears, and `poly_eval()` evaluates the
polynomial with a function that returns the observed derivative of a
tuple.

## See also

[`dexpected_by_method()`](https://statmodels7.github.io/distributions7/reference/dexpected_by_method.md)
