# The Scalar Code of a Ready-Made Transformer

Identifies a transformer built by one of the twelve ready-made
constructors and returns its code in the compiled registry and its
parameters. A candidate is rebuilt from the parameters found in the
transformer's closures, and the transformer is recognized when every one
of its functions has the candidate's code and, where it is a closure of
the constructor, the transformer's own environment; any other
transformer has no code.

## Usage

``` r
transformer_scalar_code(tr)
```

## Arguments

- tr:

  A
  [`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
  object.

## Value

`NULL`, or a list with `code`, a single integer, and `par`, a named list
of the transformer's parameters.
