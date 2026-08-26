# Higher Derivatives of the Elastic Net's Log Mills Ratio

Adds \\G''\\ and \\G'''\\ to the list of pieces, obtained by
differentiating the identity \\G' = 1 + xG - G^2\\ instead of the ratio
itself: \$\$G'' = G + xG' - 2GG', \qquad G''' = 2G' + xG'' - 2(G')^{2} -
2GG''.\$\$

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
each of the length `g` has.

## Details

Written this way each order is a polynomial in \\x\\, \\G\\ and the
orders below, so the cancellation that afflicts \\G = x - 1/M\\ for
large \\x\\ is confined to \\G\\ itself, where `.enet_G()` already
switches to an asymptotic series past \\\|x\| = 10^3\\. Nothing above
first order introduces a new cancellation.

## Notation

\\x = a/\sqrt c\\, \\M\\ the Mills ratio and \\G = \mathrm{d}\log
M/\mathrm{d}x\\.

## See also

[`enet_distrib()`](https://statmodels7.github.io/distributions7/reference/enet_distrib.md)
for the family and
[`distrib_deriv3.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.EnetDistrib.md)
for the order these serve.

## Examples

``` r
p <- distributions7:::.enet_g_higher(
  distributions7:::.enet_parts(list(mu = 0, lambda = 2, alpha = 0.5)))
c(x = p$x, G = p$g, dG = p$dg, d2G = p$d2g, d3G = p$d3g)
#>           x           G          dG         d2G         d3G 
#>  1.00000000 -0.52513528  0.19909767 -0.11693120  0.07917498 

# G'' against a central difference of the identity for G'.
gp <- function(x) { g <- distributions7:::.enet_G(x); 1 + x * g - g^2 }
eps <- 1e-5
c(analytic = p$d2g, numeric = (gp(p$x + eps) - gp(p$x - eps)) / (2 * eps))
#>   analytic    numeric 
#> -0.1169312 -0.1169312 
```
