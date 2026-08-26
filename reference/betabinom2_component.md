# One Component of a Beta-Binomial Derivative in the Shapes

Returns the derivative component of a given order and multi-index for
the beta-binomial log-mass in its two shapes. The log-mass is a sum of
log-gamma terms, so its derivative of order \\k+1\\ is the same sum with
\\\psi^{(k)}\\ in place of \\\log\Gamma\\. All four orders come from
this one routine.

## Usage

``` r
betabinom2_component(y, a, b, n, k, i, j)
```

## Arguments

- y:

  A numeric vector of counts.

- a, b:

  The two shapes, each a numeric vector of length 1 or of the length of
  `y`, strictly positive.

- n:

  The size, a single positive integer.

- k:

  The polygamma order, one less than the derivative order: 0 for the
  score, 1 for the Hessian, and so on.

- i:

  The number of \\a\\ indices in the component.

- j:

  The number of \\b\\ indices. `i + j` is the derivative order.

## Value

A numeric vector of the recycled length of the inputs. When both `i` and
`j` are positive the value is constant along it, a mixed component not
depending on `y`.

## Details

With \\a\\ and \\b\\ the shapes and \\n\\ the size, the log-mass is
\$\$\log\Gamma(y+a) + \log\Gamma(n-y+b) - \log\Gamma(n+a+b) -
\log\Gamma(a) - \log\Gamma(b) + \log\Gamma(a+b)\$\$ up to a term free of
the parameters. A derivative in \\a\\ alone differentiates the first,
third, fourth and sixth terms; one in \\b\\ alone the second, third,
fifth and sixth; a **mixed** one only the two terms carrying \\a+b\\,
which is why a mixed component of any order is \\-\psi^{(k)}(n+a+b) +
\psi^{(k)}(a+b)\\ and free of the data.

## See also

[`betabinom2_derivs()`](https://statmodels7.github.io/distributions7/reference/betabinom2_derivs.md),
which calls this for every multi-index of an order, and
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for the family.
