# Distinct Entries of a Symmetric Matrix, and Their Labels

Returns the row and column indices of the lower triangle of a \\p \times
p\\ matrix, diagonal included, together with a label per entry. A
symmetric matrix has `p * (p + 1) / 2` distinct entries, and reporting
all \\p^2\\ would give every off-diagonal quantity twice.

## Usage

``` r
mv_entry_index(p, prefix)
```

## Arguments

- p:

  The side of the matrix, a single positive whole number.

- prefix:

  The label prefix, a single string. Entry \\(i, j)\\ is named
  `prefix_vi_vj`, so `prefix = "cor"` gives `cor_v2_v1`.

## Value

A named list with `i` and `j`, integer vectors of length
`p * (p + 1) / 2` holding the row and column of each entry in
column-major order over the lower triangle, and `name`, the character
vector of labels.

## See also

[`mv_sd_cor()`](https://statmodels7.github.io/distributions7/reference/mv_sd_cor.md)
and
[`mv_derived.multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/mv_derived.multivariate_distrib.md),
which label their quantities with this.

## Examples

``` r
distributions7:::mv_entry_index(3, "cor")
#> $i
#> [1] 1 2 3 2 3 3
#> 
#> $j
#> [1] 1 1 1 2 2 3
#> 
#> $name
#> [1] "cor_v1_v1" "cor_v2_v1" "cor_v3_v1" "cor_v2_v2" "cor_v3_v2" "cor_v3_v3"
#> 

# The count is p(p + 1) / 2, the distinct entries of a symmetric matrix.
c(got = length(distributions7:::mv_entry_index(4, "sigma")$i),
  expected = 4 * 5 / 2)
#>      got expected 
#>       10       10 

# Every index pair is on or below the diagonal.
idx <- distributions7:::mv_entry_index(4, "sigma")
all(idx$i >= idx$j)
#> [1] TRUE
```
