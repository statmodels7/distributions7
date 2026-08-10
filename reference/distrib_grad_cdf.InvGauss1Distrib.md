# Inverse-Gaussian Log-CDF Derivatives

Closed form at every order. The distribution function is \\\Phi(a) +
e^{c}\Phi(b)\\ with \\a = (q/\mu - 1)/\sqrt{\phi q}\\, \\b = -(q/\mu +
1)/\sqrt{\phi q}\\ and \\c = 2/(\phi\mu)\\. Each of the three is a
product of a function of the mean and a function of the dispersion, so
its mixed partial derivatives are products of one-variable ones, and the
two terms are then a Leibniz split between the weight and the tail. The
weight and the tail are combined on the log scale: \\e^{c}\\ overflows
exactly where \\\Phi(b)\\ underflows.

## Arguments

- distrib:

  An `InvGauss1Distrib` object.

- q:

  A numeric vector of quantiles.

- theta:

  A list containing `mu` and `phi`.

- lower.tail:

  Logical; if `TRUE` (default), the lower tail.

- log:

  Logical; if `TRUE` (default), derivatives of the log probability.

- ...:

  Unused.

## Value

A named list, one vector per component.

## See also

[`invgauss1_distrib`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md)
