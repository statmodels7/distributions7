# The Pieces a von Mises Derivative in the Resultant Length Needs

Evaluates the concentration \\\kappa = A^{-1}(\rho)\\, its four
derivatives in \\\rho\\, and the derivatives of \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\ at that concentration, once per call, so that
a density or a derivative method shares them.

## Usage

``` r
vm2_parts(theta)
```

## Arguments

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1 or of a common length. `rho` must lie in \\(0, 1)\\; only it
  is read here, `mu` not entering the map.

## Value

A named list with `kappa`, the concentration; `kd`, the result of
[`numericals7::bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.html),
carrying `kappa` and its derivatives `d1` to `d4` in \\\rho\\; and `ad`,
the result of
[`numericals7::bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.html)
at that concentration, carrying `A` and its derivatives `d1` to `d4` in
\\\kappa\\.

## Details

The inverse of \\A\\ has no closed form.
[`numericals7::bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.html)
obtains \\\kappa\\ by root finding on \\\log\kappa\\ and differentiates
it by the inverse function rule, so \\\mathrm{d}\kappa/\mathrm{d}\rho =
1/A'(\kappa)\\ and the higher derivatives follow from \\A'\\ to
\\A^{(4)}\\. Those come from the Bessel recurrences and from the same
two evaluations \\A\\ already needs, so no further Bessel call is made
at any order.

The map is very steep near \\\rho = 1\\: measured,
\\\mathrm{d}\kappa/\mathrm{d}\rho\\ is 2.00 at \\\rho = 0.01\\, 3.14 at
0.5, 199 at 0.95 and \\5.0\times10^{5}\\ at 0.999.

## See also

[`distrib_gradient.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.VonMises2Distrib.md)
for the first consumer,
[`numericals7::bessel_i_ratio_inverse()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_inverse.html)
for the root finding, and
[`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md)
for the family.
