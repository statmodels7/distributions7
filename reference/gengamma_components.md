# Derivative Components of the Generalized Gamma

Returns the components of \\\partial^{\alpha+\beta+\gamma}\ell /
\partial a^\alpha \partial d^\beta \partial p^\gamma\\ at any order from
one to four, assembled term by term from the five pieces of the
log-density.

## Usage

``` r
gengamma_components(y, theta, order)
```

## Arguments

- y:

  A numeric vector of positive observations.

- theta:

  A named list with components `a`, `d` and `p`, each a numeric vector
  of length 1 or of the length of `y`, all strictly positive. Shorter
  components are recycled to the common length.

- order:

  The derivative order, an integer from 1 to 4.

## Value

A named list of component vectors, one per distinct multi-index of the
given order and keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them: three at order 1, six at order 2, ten at order 3 and fifteen
at order 4. Each has the recycled length of the inputs.

## Why the assembly is short

With \\L = \log(y/a)\\ the log-density splits into \$\$\ell = \log p -
d\log a - \log\Gamma(d/p) + (d-1)\log y - e^{pL},\$\$ and each piece is
either elementary or a univariate function composed with a
**two-variable** inner map, which the written-out template of
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
covers.

What keeps the sum from growing is that the two compositions do not
share a variable pair: \\-\log\Gamma(d/p)\\ involves \\(d, p)\\ and
\\-e^{pL}\\ involves \\(a, p)\\. Every component of the three-variable
derivative is therefore one term of one composition plus the elementary
pieces, and no genuinely three-variable expansion is ever formed. A
component naming both \\a\\ and \\d\\ comes from the elementary \\-d\log
a\\ alone.

The assembly is checked at order two against the compiled Hessian, which
was written independently, so the orders that cannot be checked against
a hand-written form rest on the orders that can.

## See also

[`distrib_deriv3.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
and
[`distrib_deriv4.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GenGamma1Distrib.md),
which call this;
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
for the two-variable composition template; and
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
for the family.
