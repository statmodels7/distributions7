# Invert the Hessian Component Names

Maps each name produced by
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
back to the pair of parameters it differentiates with respect to.

## Usage

``` r
hess_pairs(params)
```

## Arguments

- params:

  A character vector of parameter names.

## Value

A named list, parallel to `hess_names(params)`, each element a character
pair.

## Details

The wrappers need to go from `"mu_sigma"` back to `c("mu", "sigma")` in
order to combine the parent's score with its Hessian. Splitting the
string on `"_"` is the obvious way and it is wrong: a parameter whose
own name contains an underscore (`"log_scale"`) makes
`"log_scale_log_scale"` split into four pieces, and taking the first and
the last silently yields `c("log", "scale")`. Building the map from the
parameter vector cannot be fooled.

[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
is the same idea for orders above two.

## See also

[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`deriv_indices()`](https://statmodels7.github.io/distributions7/reference/deriv_indices.md)
