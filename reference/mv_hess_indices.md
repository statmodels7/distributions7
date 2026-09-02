# Index Pairs Behind the Hessian Keys of a Multivariate Distribution

Returns, for each key of
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
the two positions in `distrib@params` that the key names, in the same
order the keys come in. A multivariate method branches on whether a
component belongs to the mean block, the matrix block or the mixed
block, which is a question about position, so positions are what it asks
for.

## Usage

``` r
mv_hess_indices(distrib)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object, or any `distrib` whose `params` are the parameter names.

## Value

A list of integer vectors of length 2, one per unordered pair, as long
as `hess_names(distrib@params)` and in that order. Each pair is `(a, b)`
with `a` and `b` positions in `distrib@params`.

## Details

[`hess_pairs()`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md)
answers the same question with parameter names, which suits a univariate
method looking a closed-form component up in a table. Here the names are
matched back to positions, so the two orderings come from one
enumeration and cannot drift apart.

## See also

[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
for the keys,
[`hess_pairs()`](https://statmodels7.github.io/distributions7/reference/hess_pairs.md)
for the same enumeration as names, and
[`distrib_hessian.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MvGaussianDistrib.md)
for the consumer.

## Examples

``` r
d <- mvgaussian1_distrib(2)
idx <- distributions7:::mv_hess_indices(d)
length(idx) == length(hess_names(d@params))
#> [1] TRUE

# Key by key, the positions name the parameters the key is built from.
head(data.frame(key = hess_names(d@params),
                a = vapply(idx, `[`, integer(1), 1L),
                b = vapply(idx, `[`, integer(1), 2L)), 6)
#>                                                 key a b
#> mu1_mu1                                     mu1_mu1 1 1
#> mu2_mu2                                     mu2_mu2 2 2
#> sigma_log_L1_sigma_log_L1 sigma_log_L1_sigma_log_L1 3 3
#> sigma_log_L2_sigma_log_L2 sigma_log_L2_sigma_log_L2 4 4
#> sigma_L2.1_sigma_L2.1         sigma_L2.1_sigma_L2.1 5 5
#> mu1_mu2                                     mu1_mu2 1 2

# Which is the same enumeration hess_pairs() gives as names.
all.equal(lapply(idx, function(p) d@params[p]),
          distributions7:::hess_pairs(d@params),
          check.attributes = FALSE)
#> [1] TRUE
```
