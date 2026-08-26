# The Two Preimages of a Folded Point

Computes the parent's density at \\+x\\ and at \\-x\\, their sum, and
the weight \\w = f(x)/L(x)\\ of the positive preimage. Every method of
[FoldedDistrib](https://statmodels7.github.io/distributions7/reference/FoldedDistrib.md)
needs the same four quantities, so they are formed once per call: the
density is \\L\\, the score is \\w s(x) + (1-w) s(-x)\\, and every
higher-order ratio is a \\w\\-weighted average of the parent's own.

## Usage

``` r
fold_parts(parent, x, theta)
```

## Arguments

- parent:

  The wrapped `continuous_distrib` object.

- x:

  A numeric vector of points at which to evaluate. Negative values are
  not screened out here; the calling method zeroes them.

- theta:

  A named list of the parent's parameters.

## Value

A named list of four numeric vectors of the recycled length: `fp` and
`fm`, the parent's density at \\+x\\ and \\-x\\; `L`, their sum; and
`w`, the ratio `fp / L`.

## Details

At \\x = 0\\ the two preimages coincide, so \\w = 1/2\\ whatever the
parameters are and the folded density is exactly twice the parent's. Far
out in the tail of a parent centered above zero, \\w\\ approaches one
and the fold becomes invisible.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one, \\w\\ the weight of
the positive preimage and \\s\\ the parent's score.

## See also

[`fold_ratio()`](https://statmodels7.github.io/distributions7/reference/fold_ratio.md)
for the higher-order quantities built on these, and
[`distrib_pdf.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.FoldedDistrib.md)
for the first consumer.

## Examples

``` r
g <- gaussian1_distrib()
theta <- list(mu = 0.5, sigma = 1.2)
p <- distributions7:::fold_parts(g, c(0, 1, 3), theta)
str(p)
#> List of 4
#>  $ fp: num [1:3] 0.305 0.305 0.038
#>  $ fm: num [1:3] 0.30481 0.15221 0.00473
#>  $ L : num [1:3] 0.6096 0.457 0.0427
#>  $ w : num [1:3] 0.5 0.667 0.889

# L is the sum of the two preimages, which is the folded density.
all.equal(p$L, dnorm(c(0, 1, 3), 0.5, 1.2) + dnorm(-c(0, 1, 3), 0.5, 1.2))
#> [1] TRUE

# At zero the two preimages coincide, so the weight is exactly one half.
p$w[1]
#> [1] 0.5

# And far into the tail the fold becomes invisible: w approaches one.
distributions7:::fold_parts(g, c(3, 6, 12), theta)$w
#> [1] 0.8892727 0.9847328 0.9997597
```
