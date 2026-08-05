# Higher Derivatives of the Beta in Its Shapes

Orders two to four, which are free of the data: each is a difference of
polygamma functions at \\\alpha\\, \\\beta\\ and \\\alpha+\beta\\.

## Usage

``` r
beta2_higher(theta, n, order)
```

## Arguments

- theta:

  A list with `alpha` and `beta`.

- n:

  The number of observations to recycle to.

- order:

  The derivative order, 2, 3 or 4.

## Value

A named list of component vectors.

## Details

Being free of the data, the observed and the expected derivatives
coincide at every order beyond the first, which is why this family needs
no expectation anywhere.

## See also

[`beta2_distrib`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md)
