# Split a Wrapper's Parameters From Its Parent's

Separates the full `theta` into the parent distribution's parameters and
the single mixture parameter the wrapper adds.

## Usage

``` r
split_mix_theta(distrib, theta)
```

## Arguments

- distrib:

  A zero-inflated or zero-adjusted distribution object.

- theta:

  A named list of parameters, already aligned.

## Value

A list with `orig`, the parent's parameters, and `mix`, the wrapper's
own.

## Details

Both zero wrappers append their parameter last, so the split is
positional and does not depend on what that parameter is called – `zi`
for inflation, `pi` for adjustment.
