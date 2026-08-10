# The Skew Normal's Distribution Function in Standard Coordinates

Returns \\\partial^{i}\_{z}\partial^{j}\_{\alpha}G\\ for \\G(z, \alpha)
= \Phi(z) - 2T(z, \alpha)\\, over the pairs with \\1 \le i + j \le\\
`order`, as a list indexed by `i + 1` then `j + 1`.

## Usage

``` r
sn_cdf_std_derivs(z, al, order)
```

## Arguments

- z:

  The standardized quantile.

- al:

  The shape.

- order:

  The highest total order, 1 to 4.

## Value

A nested list, `[[i + 1]][[j + 1]]`.

## Details

The first derivatives are \\G\_{z} = 2\varphi(z)\Phi(u)\\ and
\\G\_{\alpha} = -2\varphi(z)\varphi(u)R\\ with \\u = \alpha z\\ and \\R
= (1+\alpha^{2})^{-1}\\; everything above them is one Leibniz rule over
those factors, with \\\Phi(u)\\ and \\\varphi(u)\\ differentiated by Faa
di Bruno over \\u\\, which is bilinear and so has only three non-zero
partial derivatives.

## See also

[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
