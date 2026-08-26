# The Block Ratios of a Folded Density

Returns a memoized function giving \\d^B L / L\\ for any block \\B\\ of
parameter names. That is what
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
consumes to assemble a derivative of \\\log L\\ of any order. The ratio
is the parent's complete Bell polynomial evaluated at each preimage and
averaged with the weight \\w\\, \$\$\frac{d^B L}{L} = w \frac{d^B
f(x)}{f(x)} + (1-w) \frac{d^B f(-x)}{f(-x)}.\$\$

## Usage

``` r
fold_ratio(parent, x, theta, order, params, w)
```

## Arguments

- parent:

  The wrapped `continuous_distrib` object.

- x:

  A numeric vector of points.

- theta:

  A named list of the parent's parameters.

- order:

  The highest derivative order the caller will ask for, a single whole
  number from 1 to 4. Blocks longer than this are not prepared.

- params:

  The parameter names, as `distrib@params`.

- w:

  The weight of the positive preimage, from
  [`fold_parts()`](https://statmodels7.github.io/distributions7/reference/fold_parts.md).

## Value

A function of one character vector, the block, returning a numeric
vector of the recycled length of `x`.

## Details

The memoization matters because a partition sum asks for the same block
many times: at order four the same singleton block appears in most of
the partitions. Each distinct block costs two evaluations of the
parent's derivatives, one at \\+x\\ and one at \\-x\\, and no more.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one, \\B\\ a multiset of
parameter names, \\d^B\\ the corresponding partial derivative and \\w\\
the weight of the positive preimage.

## See also

[`fold_parts()`](https://statmodels7.github.io/distributions7/reference/fold_parts.md)
for the weight,
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
for the partition sum this feeds, and
[`bell_f_ratio()`](https://statmodels7.github.io/distributions7/reference/bell_f_ratio.md)
for the parent's own Bell quantity.

## Examples

``` r
g <- gaussian1_distrib()
theta <- list(mu = 0.5, sigma = 1.2)
x <- c(0.5, 1, 3)
w <- distributions7:::fold_parts(g, x, theta)$w
r <- distributions7:::fold_ratio(g, x, theta, 2L, c("mu", "sigma"), w)

# A singleton block is the weighted average of the parent's own score.
sp <- distrib_gradient(g, x, theta)$mu
sm <- distrib_gradient(g, -x, theta)$mu
all.equal(r("mu"), w * sp + (1 - w) * sm)
#> [1] TRUE

# Which at first order IS the folded score, log L having no other term.
all.equal(r("mu"), distrib_gradient(folded(g), x, theta)$mu)
#> [1] TRUE

# A two-element block is the second-order ratio, not the second derivative
# of the logarithm: the two differ by the square of the first.
c(ratio = r(c("mu", "mu"))[2],
  log_second = distrib_hessian(folded(g), x, theta)$mu_mu[2],
  difference = r(c("mu", "mu"))[2] - r("mu")[2]^2)
#>      ratio log_second difference 
#> -0.2526571 -0.2659605 -0.2659605 
```
