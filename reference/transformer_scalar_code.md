# The Scalar Code of a Ready-Made Transformer

Identifies a transformer built by one of the twelve ready-made
constructors and returns its code in the compiled registry and its
parameters. A candidate is rebuilt from the parameters found in the
transformer's closures, and the transformer is recognized when every one
of its eight functions has the same source text as the candidate's, with
the same `decreasing` flag, and, where the candidate's function is a
closure created by the constructor, the same environment as the
transformer's own function; any other transformer has no code.

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
