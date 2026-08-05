# The Concentration a Mean Resultant Length Implies

\\\kappa = A^{-1}(\rho)\\, by root finding, together with the four
derivatives of the inverse.

## Usage

``` r
vm_kappa_of_rho(rho)
```

## Arguments

- rho:

  The mean resultant length, in \\(0, 1)\\.

## Value

A named list with `kappa` and its four derivatives in `rho`.

## Details

\\A\\ has no elementary inverse, so \\\kappa\\ is found by bisection on
\\\log\kappa\\, where the function is well conditioned over the whole
range. The derivatives then come from the inverse function rule, which
needs no further root finding: \$\$\kappa' = \dfrac{1}{A'}, \qquad
\kappa'' = -\dfrac{A''}{(A')^3}, \qquad \kappa''' = \dfrac{3(A'')^2 -
A'A'''}{(A')^5},\$\$ and the fourth in the same pattern. \\A'\\ is the
variance of \\\cos(Y-\mu)\\ and therefore strictly positive, so none of
these divides by zero in the interior.

## See also

[`vonmises2_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
