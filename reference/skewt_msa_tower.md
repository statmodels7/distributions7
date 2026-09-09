# The Skew t Tower in the Location, Scale and Shape

Builds the table of \\\partial^i_z\\\partial^c\_\alpha \ell\\ that the
third and fourth derivatives of a skew t read, for every pair with \\c +
i \le 4\\. Together with
[`skewt_msa_component()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_component.md)
it gives in closed form every derivative of the log-density that does
not involve \\\nu\\.

## Usage

``` r
skewt_msa_tower(y, mu, sigma, alpha, nu)
```

## Arguments

- y:

  A numeric vector of observations.

- mu, sigma, alpha, nu:

  The four parameters, each of length 1 or of the length of `y`.

## Value

A named list. `z` is the standardized residual; the remaining fifteen
elements are named `"c_i"` and hold \\\partial_z^i \Phi_c(z)\\ for every
\\c + i \le 4\\, with \\\Phi_0 = \ell\\ up to the terms free of \\z\\,
so that `"0_0"` is `NA_real_` and is never read.

## Why the block closes

With \\z = (y-\mu)/\sigma\\ and \\u(z) = z\sqrt{(\nu+1)/(\nu+z^2)}\\,
the log-density is \$\$\ell = \log 2 - \log\sigma + g(z) +
\Lambda(\alpha\\u(z)),\$\$ where \\g = \log t\_\nu\\ and \\\Lambda =
\log T\_{\nu+1}\\. Neither \\\mu\\ nor \\\sigma\\ appears anywhere but
inside \\z\\, and \\\alpha\\ nowhere but inside \\w = \alpha u\\, so two
observations settle the whole block.

The shape enters the argument of \\\Lambda\\ linearly, so
differentiating in it never leaves the table: \$\$\Phi_c(z) \\=\\
\dfrac{\partial^c \ell}{\partial\alpha^c} \\=\\
u(z)^c\\\Lambda^{(c)}(w), \qquad c \ge 1 .\$\$

The location and the scale then act on a function of \\z\\ alone, and
for any such \\F\\, \$\$\dfrac{\partial^{a+b}
F}{\partial\mu^a\\\partial\sigma^b} =
\dfrac{(-1)^{a+b}}{\sigma^{a+b}}\\P\_{a,b}(z), \qquad P\_{a,b}(z) =
\sum\_{i=0}^{b}\binom{b}{i} \left\[\prod\_{k=i}^{b-1}(a+k)\right\]
z^i\\F^{(a+i)}(z),\$\$ which follows by induction from \\P\_{a,0} =
F^{(a)}\\ and \\P\_{a,b+1} = (a+b)\\P\_{a,b} + z\\P\_{a,b}'\\. The
explicit \\-\log\sigma\\ contributes only to the pure-\\\sigma\\
components, where it adds \\(-1)^b (b-1)!\\\sigma^{-b}\\.

## Nothing in it needs a difference

\\u\\ and \\g\\ are rational in \\z\\ up to one square root, so their
four derivatives are written out. \\\Lambda^{(k)}\\ follows from the
Riccati recursion \\Q' = Q(G-Q)\\ for \\Q = t\_{\nu+1}/T\_{\nu+1}\\ and
\\G = \partial_w \log t\_{\nu+1}\\, differentiated twice more, exactly
as a skew normal's inverse Mills ratio is handled. The two chain rules
that remain, \\(u^c)^{(j)}\\ and \\\partial_z^p \Lambda^{(c)}(\alpha
u)\\, are Faa di Bruno sums over the partial Bell polynomials of \\u\\.

What is left outside the table is \\\nu\\, which enters the degrees of
freedom of \\T\_{\nu+1}\\ and therefore carries the obstruction
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md)
records.

## See also

[`skewt_msa_component()`](https://statmodels7.github.io/distributions7/reference/skewt_msa_component.md),
which reads one derivative off the table, and
[`skewt_pieces()`](https://statmodels7.github.io/distributions7/reference/skewt_pieces.md),
whose `a`, `e` and `q` are this function's \\g'\\, \\u'\\ and \\Q\\.

## Examples

``` r
tw <- distributions7:::skewt_msa_tower(c(-0.4, 1.2), 0, 1, 0.7, 8)
names(tw)
#>  [1] "z"   "0_0" "0_1" "0_2" "0_3" "0_4" "1_0" "1_1" "1_2" "1_3" "2_0" "2_1"
#> [13] "2_2" "3_0" "3_1" "4_0"
```
