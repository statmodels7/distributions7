# Index Tuples of a Given Width Over a Number of Variables

Enumerates the index multisets of exactly `order` indices drawn from
`1:d`, in the order
[`parameters7::param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.html)
uses. It takes a COUNT rather than a parametrization, so a composite
index set formed by appending a coordinate can be enumerated without
building an object of that shape. The Dirichlet's set is the one that
needs this: its simplex free values followed by its concentration. The
work is done by
[`numericals7::tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.html),
the toolkit's one copy of the enumeration.

## Usage

``` r
tuple_indices_upto(d, order)
```

## Arguments

- d:

  The number of variables, a single positive whole number.

- order:

  The tuple width, a single whole number from 1 to 4.

## Value

A list of integer vectors, each of length `order`, sorted within itself.
Its length is `choose(d + order - 1, order)`.

## See also

[`numericals7::tuple_indices()`](https://statmodels7.github.io/numericals7/reference/tuple_indices.html)
for the enumeration itself,
[`dirichlet_map_tensors()`](https://statmodels7.github.io/distributions7/reference/dirichlet_map_tensors.md)
for the consumer, and
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
for the same enumeration rendered as component names.

## Examples

``` r
# Two indices over three variables: the three repeats then the three pairs.
distributions7:::tuple_indices_upto(3, 2)
#> [[1]]
#> [1] 1 1
#> 
#> [[2]]
#> [1] 2 2
#> 
#> [[3]]
#> [1] 3 3
#> 
#> [[4]]
#> [1] 1 2
#> 
#> [[5]]
#> [1] 1 3
#> 
#> [[6]]
#> [1] 2 3
#> 

# The count is choose(d + order - 1, order).
c(got = length(distributions7:::tuple_indices_upto(4, 3)),
  expected = choose(4 + 3 - 1, 3))
#>      got expected 
#>       20       20 

# Every tuple is sorted, which is what makes a comma-joined key canonical.
all(vapply(distributions7:::tuple_indices_upto(3, 3),
           function(t) !is.unsorted(t), TRUE))
#> [1] TRUE
```
