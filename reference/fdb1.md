# Univariate Composition, One Order

The derivative of order `order` of \\h(u(x))\\ for scalar chains, Faa di
Bruno written out: \\(h \circ u)'' = h''u_1^2 + h'u_2\\ and so on to
order four. Each order is its own formula and forms nothing above it.

## Usage

``` r
fdb1(h, u, order)
```

## Arguments

- h:

  A list of the outer derivatives `h1` to at least `h_order` at the
  inner value.

- u:

  A list of the inner derivatives `u1` to at least `u_order`.

- order:

  An integer from 1 to 4.

## Value

A numeric vector, the derivative of order `order` of the composition.

## See also

[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
for a bivariate inner function;
[`reparam_map_derivs()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md),
which consumes it.
