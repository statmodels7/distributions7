# Partial Derivatives of a Product of a Location Term and a Scale Term

Builds the `wderiv` or `xderiv` callback
[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
wants, for a quantity that factorizes as \\u(\theta_1)\\v(\theta_2)\\. A
mixed partial of such a product is the product of two one-variable
derivatives, so a block of parameter names is answered by counting how
many of each it holds.

## Usage

``` r
separable_deriv(nm, uderiv, vderiv)
```

## Arguments

- nm:

  A character vector of length 2, the two parameter names, the one
  `uderiv` differentiates first.

- uderiv:

  A function of a non-negative whole number \\j\\ returning \\\partial^j
  u/\partial\theta_1^j\\.

- vderiv:

  The same for \\v\\ in the second parameter.

## Value

A function of a block, a character vector of parameter names, that
returns the corresponding mixed partial. Names outside `nm` are not
counted, so a block made only of such names gives the order-zero factor
for both, which is the value itself.

## Details

The inverse Gaussian, which uses this function, is separable in this
way. The inverse Gaussian's \\a\\, \\b\\ and \\c\\ are each a function
of the mean times a function of the dispersion; the elastic net's \\s\\
and \\x\\ are each a function of \\\lambda\\ times a function of
\\\alpha\\. Each mixed partial is therefore a product of two
one-variable derivatives, which is why four orders are cheap here.

## See also

[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md),
whose `wderiv` and `xderiv` this builds;
[`dpow_affine()`](https://statmodels7.github.io/distributions7/reference/dpow_affine.md)
for the other shape.
