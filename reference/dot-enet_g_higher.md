# Higher Derivatives of the Elastic Net's Log Mills Ratio

\\G''\\ and \\G'''\\, obtained by differentiating the identity \\G' =
1 + xG - G^2\\ rather than the ratio itself.

## Usage

``` r
.enet_g_higher(p)
```

## Arguments

- p:

  The value of `.enet_parts`.

## Value

`p` with `d2g` and `d3g` added.

## Details

\$\$G'' = G + xG' - 2GG', \qquad G''' = 2G' + xG'' - 2(G')^{2} -
2GG''.\$\$ Written this way each order is a polynomial in \\x\\, \\G\\
and the orders below, so the cancellation that afflicts \\G = x - 1/M\\
for large \\x\\ is confined to \\G\\ itself, where `.enet_G` already
switches to an asymptotic series.
