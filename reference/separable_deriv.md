# Partial Derivatives of a Product of a Location Term and a Scale Term

Builds the evaluator for a quantity \\U(\mu)V(\phi)\\, whose mixed
partial is the product of the two one-variable derivatives.

## Usage

``` r
separable_deriv(nm, uderiv, vderiv)
```

## Arguments

- nm:

  The two parameter names, in order.

- uderiv:

  A function of the order returning \\U^{(j)}\\.

- vderiv:

  A function of the order returning \\V^{(k)}\\.

## Value

A function of a character vector of parameter names.
