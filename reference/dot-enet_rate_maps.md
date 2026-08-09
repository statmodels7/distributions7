# The Map From the Rates to the Elastic Net's Parameters

The partial tables of \\(\mu, a, c)\\ as functions of \\(\mu, \lambda,
\alpha)\\, with \\a = \lambda\alpha\\ and \\c = \lambda(1-\alpha)\\.

## Usage

``` r
.enet_rate_maps(p)
```

## Arguments

- p:

  The value of `.enet_parts`.

## Value

A list of keyed tables, one per rate coordinate.

## Details

The map is bilinear, so its second partials are the constants \\\pm 1\\
and its third and higher vanish. Those second partials are what carries
a low-order derivative in the rates up to a high-order one in the
parameters, which is why the third derivatives in \\(\lambda, \alpha)\\
still depend on the data although the third derivatives in the rates do
not.
