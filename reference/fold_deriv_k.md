# Build a Folded Derivative Method of a Given Order

Returns the method
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
registers for
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
or
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md).
Both orders run the same three steps, so the body is written once and
the order is closed over: form the weight with
[`fold_parts()`](https://statmodels7.github.io/distributions7/reference/fold_parts.md),
build the block ratios with
[`fold_ratio()`](https://statmodels7.github.io/distributions7/reference/fold_ratio.md),
and hand them to
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
for every component of the order.

## Usage

``` r
fold_deriv_k(order)
```

## Arguments

- order:

  The derivative order, `3L` or `4L`.

## Value

A function with the signature of
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
suitable for `S7::method(...) <- `.

## Details

The ratios are \\d^B L/L = w\\(d^B f(x)/f(x)) + (1-w)\\(d^B
f(-x)/f(-x))\\, and
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
turns them into derivatives of \\\log L\\ by the moment-to-cumulant
relation. Nothing about the fold is written out at third or fourth
order; the same two functions serve every order the parent supplies.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one, \\w\\ the weight of
the positive preimage and \\B\\ a multiset of parameter names.

## See also

[`distrib_deriv3.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.FoldedDistrib.md)
and
[`distrib_deriv4.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.FoldedDistrib.md),
the two methods it builds, and
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
for the partition sum.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)
set.seed(1)
y <- distrib_rng(d, 6, theta)

m3 <- distributions7:::fold_deriv_k(3L)

# It builds the method the class registers, so the two agree.
all.equal(m3(d, y, theta), distrib_deriv3(d, y, theta))
#> [1] TRUE

# And the order it was built at fixes the component count.
c(order3 = length(m3(d, y, theta)),
  order4 = length(distributions7:::fold_deriv_k(4L)(d, y, theta)))
#> order3 order4 
#>      4      5 
```
