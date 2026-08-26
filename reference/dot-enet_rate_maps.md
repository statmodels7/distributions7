# The Map From the Rates to the Elastic Net's Parameters

Returns the partial derivatives of \\(\mu, a, c)\\ as functions of
\\(\mu, \lambda, \alpha)\\, with \\a = \lambda\alpha\\ and \\c =
\lambda(1-\alpha)\\, in the keyed-table form
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
consumes.

## Usage

``` r
.enet_rate_maps(p)
```

## Arguments

- p:

  The value of `.enet_parts()`, read for `al` and `lam`.

## Value

A named list of three keyed tables, one per rate coordinate, with keys
`"1"`, `"2"`, `"3"` for the first partials in \\\mu\\, \\\lambda\\,
\\\alpha\\ and `"2,3"` for the one non-zero second partial.

## Details

The map is bilinear, so its second partials are the constants \\\pm 1\\
and its third and higher vanish; the table says so by leaving them out.

Those second partials are what carries a low-order derivative in the
rates up to a high-order one in the parameters. That is why the third
derivatives in \\(\lambda, \alpha)\\ still depend on the data although
the third derivatives in the rates do not.

## Notation

\\a = \lambda\alpha\\, \\c = \lambda(1-\alpha)\\, and the numeric keys
index \\(\mu, \lambda, \alpha)\\ in that order.

## See also

[`.enet_chain()`](https://statmodels7.github.io/distributions7/reference/dot-enet_chain.md),
its only caller, and
[`chain_assemble()`](https://statmodels7.github.io/distributions7/reference/chain_assemble.md)
for the partition sum that consumes it.

## Examples

``` r
p <- distributions7:::.enet_parts(list(mu = 0, lambda = 2, alpha = 0.5))
distributions7:::.enet_rate_maps(p)
#> $mu
#> $mu$`1`
#> [1] 1
#> 
#> 
#> $a
#> $a$`2`
#> [1] 0.5
#> 
#> $a$`3`
#> [1] 2
#> 
#> $a$`2,3`
#> [1] 1
#> 
#> 
#> $c
#> $c$`2`
#> [1] 0.5
#> 
#> $c$`3`
#> [1] -2
#> 
#> $c$`2,3`
#> [1] -1
#> 
#> 

# The two rates add up to lambda whatever alpha is, so their first
# partials in lambda sum to one and those in alpha cancel.
m <- distributions7:::.enet_rate_maps(p)
c(in_lambda = m$a[["2"]] + m$c[["2"]],
  in_alpha = m$a[["3"]] + m$c[["3"]])
#> in_lambda  in_alpha 
#>         1         0 
```
