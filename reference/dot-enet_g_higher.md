# Higher Derivatives of the Elastic Net's Log Mills Ratio

Adds \\G''\\ and \\G'''\\ to the list of pieces, the second and third
derivatives of \\G = \mathrm{d}\log M/\mathrm{d}x\\.

## Usage

``` r
.enet_g_higher(p)
```

## Arguments

- p:

  The value of `.enet_parts()`: a list carrying `mu`, `lam`, `al`, `a`,
  `c`, `x`, `g` and `dg`.

## Value

`p` with two components added: `d2g` for \\G''\\ and `d3g` for \\G'''\\,
and each of them has the same length as `g`.

## Details

The derivatives of \\\log M\\ are, with alternating signs, the cumulants
of the density proportional to \\e^{-xs - s^2/2}\\ on \\s \> 0\\, and a
compiled kernel computes them from the ratios of its moments, which
follow a recursion run downward for \\x \ge 1\\. The identities \\G' =
1 + xG - G^2\\ and \\G'' = G + xG' - 2GG'\\ hold, but they cancel for
large \\x\\ (\\1 + xG\\ is of order \\x^{-2}\\ while its terms are of
order one), so they are not used.

## Notation

\\x = a/\sqrt c\\, \\M\\ the Mills ratio and \\G = \mathrm{d}\log
M/\mathrm{d}x\\.

## See also

[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
for the family and
[`distrib_deriv3.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md)
for the order that these serve.

## Examples

``` r
p <- distributions7:::.enet_g_higher(
  distributions7:::.enet_parts(list(mu = 0, lambda = 2, alpha = 0.5)))
c(x = p$x, G = p$g, dG = p$dg, d2G = p$d2g, d3G = p$d3g)
#>           x           G          dG         d2G         d3G 
#>  1.00000000 -0.52513528  0.19909767 -0.11693120  0.07917498 

# G'' is compared with a central difference of G'.
gp <- function(x) {
  distributions7:::.enet_parts(list(mu = 0, lambda = 2 * x^2,
                                    alpha = 0.5))$dg
}
eps <- 1e-5
c(analytic = p$d2g, numeric = (gp(p$x + eps) - gp(p$x - eps)) / (2 * eps))
#>   analytic    numeric 
#> -0.1169312 -0.1169312 
```
