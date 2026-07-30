# Assemble One Order of a Wrapper's Derivatives

The shared skeleton: builds the named list of components of a given
order by calling `component` once per multi-index.

## Usage

``` r
assemble_deriv(distrib, order, component)
```

## Arguments

- distrib:

  The wrapper distribution.

- order:

  The derivative order.

- component:

  A function of one multi-index returning that component.

## Value

A named list of derivative component vectors.
