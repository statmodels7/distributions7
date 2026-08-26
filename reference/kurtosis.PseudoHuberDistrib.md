# Excess Kurtosis of the Pseudo-Huber Distribution

Closed form, replacing the numerical default: \$\$\gamma_2 =
3\\\frac{K_3(\sqrt{\nu})\\K_1(\sqrt{\nu})} {K_2(\sqrt{\nu})^2} - 3,\$\$
with \\K_r\\ the modified Bessel function of the second kind. It depends
on the shape alone, the location and the scale canceling in a
standardized moment, and it interpolates between the two limits the
family is built to span: about 3, the Laplace's, as \\\nu \to 0\\, and
0, the Gaussian's, as \\\nu\\ grows.

## Arguments

- x:

  A `PseudoHuberDistrib`, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or `n`. Only `nu` enters the value; the other two
  are read for their lengths and are still validated.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of the length the recycled
parameters imply. Only `nu` enters the value, so a setting that varies
the location or the scale alone repeats one number.

## Details

All three Bessel functions are evaluated exponentially scaled. Each
carries the factor \\e^{\sqrt{\nu}}\\, and the ratio is arranged so that
the factors cancel exactly, so the expression stays usable at large
\\\nu\\ where the unscaled functions underflow.

## Notation

\\\nu \> 0\\ is the shape and \\K_r\\ the modified Bessel function of
the second kind of order \\r\\.

## See also

[`variance.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.PseudoHuberDistrib.md)
for the other Bessel ratio,
[`kurtosis.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.LaplaceDistrib.md)
and
[`kurtosis.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian1Distrib.md)
for the two limits,
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md)
for the family.

## Examples

``` r
d <- pseudohuber_distrib()

# The shape carries the tail weight, from the Laplace's 3 to the Gaussian's 0.
round(kurtosis(d, list(mu = 0, sigma = 1, nu = c(1e-4, 1, 1e2, 1e4))), 4)
#> [1] 2.9987 1.8570 0.2954 0.0300

# Neither the location nor the scale enters a standardized moment.
kurtosis(d, list(mu = c(0, 5), sigma = c(1, 9), nu = 4))
#> [1] 1.218426 1.218426
```
