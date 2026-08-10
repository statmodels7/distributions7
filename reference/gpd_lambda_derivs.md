# Derivatives of log1p(u)/u

Returns \\\Lambda(u) = \log(1+u)/u\\ and its first four derivatives, one
vector per order.

## Usage

``` r
gpd_lambda_derivs(u)
```

## Arguments

- u:

  A numeric vector, greater than \\-1\\.

## Value

A list of five numeric vectors, orders 0 to 4.

## Details

Differentiating \\u\Lambda = \log(1+u)\\ gives \\u\Lambda^{(r)} +
r\Lambda^{(r-1)} = (-1)^{r-1}(r-1)!/(1+u)^{r}\\, a recursion that is
exact away from the origin and useless at it: it divides by \\u\\ and
subtracts two nearly equal quantities, and measured against the Taylor
series the fourth derivative is wrong by a factor of \\10^{39}\\ at \\u
= 10^{-14}\\, by 1.7 at \\10^{-4}\\ and by \\3\times10^{-8}\\ at
\\10^{-2}\\. The series \\\Lambda^{(r)}(u) = \sum\_{m\ge
r}(-1)^{m}\frac{m!}{(m-r)!} \frac{u^{m-r}}{m+1}\\ is used instead below
\\\lvert u\rvert = 1/2\\, where the two agree to \\10^{-16}\\, and its
truncation is set so that the \\m^{4}\\ weight at the switch point stays
under the rounding.

## See also

[`gpd_surv_pieces`](https://statmodels7.github.io/distributions7/reference/gpd_surv_pieces.md)
