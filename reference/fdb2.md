# Bivariate Composition to Fourth Order

The fourteen partials of \\h(u(x, z))\\ for a scalar outer \\h\\ and a
bivariate inner \\u\\, Faa di Bruno written out component by component.
The inner partials arrive as a named list with entries `x`, `z`, `xx`,
`xz`, `zz`, `xxx`, `xxz`, `xzz`, `zzz`, `xxxx`, `xxxz`, `xxzz`, `xzzz`,
`zzzz`; missing entries count as zero.

## Usage

``` r
fdb2(h, u)
```

## Arguments

- h:

  A list of the outer derivatives `h1` to `h4` at the inner value.

- u:

  The named list of inner partials.

## Value

A named list of the fourteen partials of the composition, keyed exactly
as `u` is: `x`, `z`, `xx`, `xz`, `zz`, `xxx`, `xxz`, `xzz`, `zzz`,
`xxxx`, `xxxz`, `xxzz`, `xzzz`, `zzzz`.

## See also

[`fdb1()`](https://statmodels7.github.io/distributions7/reference/fdb1.md)
for a univariate inner function;
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md),
which consumes these through `map_derivs`;
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
for the general partition sum this writes out.
