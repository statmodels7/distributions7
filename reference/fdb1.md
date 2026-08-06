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

A list with the four derivatives of the composition.
