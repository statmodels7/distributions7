# Set Partitions of a Multi-Index

Returns every way of splitting a multi-index into blocks. The blocks
hold the parameter names and not the positions.

## Usage

``` r
index_partitions(idx)
```

## Arguments

- idx:

  A character vector of parameter names, with repetition.

## Value

A list of partitions, each of which is a list of character-vector
blocks.

## See also

[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
