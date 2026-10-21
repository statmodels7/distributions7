# Partial Derivatives of a Product of a Location Term and a Scale Term

Builds the `wderiv` or `xderiv` callback that
[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
calls, for a quantity that factorizes as \\u(\theta_1)\\v(\theta_2)\\. A
mixed partial of such a product is the product of two one-variable
derivatives, so the partial for a block of parameter names is obtained
by counting how many times each of the two names occurs in the block.

## Usage

``` r
separable_deriv(nm, uderiv, vderiv)
```

## Arguments

- nm:

  A character vector of length 2 holding the two parameter names.
  `uderiv` differentiates the first of them.

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
one-variable derivatives, so four orders are cheap here.

## See also

[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md),
whose `wderiv` and `xderiv` this builds;
[`dpow_affine()`](https://statmodels7.github.io/distributions7/reference/dpow_affine.md)
for the other shape.
