# Higher Derivatives of the Mean Resultant Length

\\A''\\, \\A'''\\ and \\A''''\\ at \\\kappa\\, obtained by
differentiating \\A' = 1 - A/\kappa - A^2\\ repeatedly.

## Usage

``` r
vm_A_derivs(kappa)
```

## Arguments

- kappa:

  The concentration.

## Value

A named list with `A` and its four derivatives.

## Details

Each order is written in the orders below it, so the whole table costs
the two Bessel functions
[`vm_A`](https://statmodels7.github.io/distributions7/reference/vm_A.md)
already evaluates and nothing more: \$\$A'' = -\dfrac{A'}{\kappa} +
\dfrac{A}{\kappa^2} - 2AA'\$\$ and so on. The alternative, evaluating a
Bessel function of higher order for each derivative, costs more and is
less accurate at large \\\kappa\\, where the functions themselves
overflow and only their ratio does not.

## See also

[`vonmises2_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md),
[`vm_A`](https://statmodels7.github.io/distributions7/reference/vm_A.md)
