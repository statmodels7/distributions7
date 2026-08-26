# Univariate Composition to Fourth Order

The derivatives of \\h(u(x))\\ for scalar chains, Faa di Bruno written
out: \\(h \circ u)'' = h''u_1^2 + h'u_2\\ and so on to order four.

## Usage

``` r
fdb1(h, u)
```

## Arguments

- h:

  A list of the outer derivatives `h1` to `h4` at the inner value.

- u:

  A list of the inner derivatives `u1` to `u4`.

## Value

A named list with components `h1` to `h4`, the four derivatives of the
composition, each the length of the vectors passed in.

## See also

[`fdb2()`](https://statmodels7.github.io/distributions7/reference/fdb2.md)
for a bivariate inner function;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
which consumes these through `map_derivs`;
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
for the general partition sum this writes out.
