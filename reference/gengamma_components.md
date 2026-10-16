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

## Structure of the assembly

With \\L = \log(y/a)\\ the log-density splits into \$\$\ell = \log p -
d\log a - \log\Gamma(d/p) + (d-1)\log y - e^{pL},\$\$ and each piece is
either elementary or a univariate function composed with a two-variable
inner map, which the written-out template of
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
covers.

The sum does not grow, because the two compositions do not share a
variable pair: \\-\log\Gamma(d/p)\\ involves \\(d, p)\\ and \\-e^{pL}\\
involves \\(a, p)\\. Every component of the three-variable derivative is
therefore one term of one composition plus the elementary pieces, and
the assembly never forms a genuinely three-variable expansion. A
component that names both \\a\\ and \\d\\ comes only from the elementary
term \\-d\log a\\.

The tests compare the assembly at orders one to four with the compiled
kernels of
[`distrib_gradient.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.GenGamma1Distrib.md),
[`distrib_hessian.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GenGamma1Distrib.md)
and the orders three and four, which were written independently. Below
\\k = d/p\\ of about 0.01 the assembly loses accuracy, since it
evaluates \\\psi^{(n)}(k)\\ near its pole.

## See also

[`distrib_deriv3.GenGamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GenGamma1Distrib.md)
for the compiled kernels against which this assembly is checked;
[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
for the two-variable composition template; and
[`gengamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gengamma1_distrib.md)
for the family.
