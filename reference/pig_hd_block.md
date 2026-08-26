# Log-Likelihood Derivatives of the Poisson-Inverse Gaussian

Evaluates the log-likelihood and its fourteen partial derivatives to
fourth order at once, through one compiled kernel, and returns the block
of components an order asks for. Both parametrizations use it, with a
different `kernel` argument; every derivative method of both families is
one call of this with a different `cols`.

## Usage

``` r
pig_hd_block(y, theta, cols, kernel, threads = 1L)
```

## Arguments

- y:

  A numeric vector of observations. Recycled with `theta` to the common
  length.

- theta:

  The parameter list, already aligned and read positionally: the mean
  first and the dispersion or \\\alpha\\ second.

- cols:

  A named character vector selecting the kernel columns wanted. The
  names become the names of the result and the values are among `"l"`,
  `"d10"`, `"d01"`, `"d20"`, `"d11"`, `"d02"`, `"d30"`, `"d21"`,
  `"d12"`, `"d03"`, `"d40"`, `"d31"`, `"d22"`, `"d13"`, `"d04"`, where
  `dij` is \\\partial^{i+j}\ell/\partial\theta_1^i\partial\theta_2^j\\.

- kernel:

  The compiled kernel: `pig1_hd_cpp` for the mean-dispersion
  parametrization or `pig2_hd_cpp` for the orthogonal one.

- threads:

  A single positive integer, how many threads that kernel may use.
  Defaults to `1L`. The result does not depend on the count.

## Value

A named list of numeric vectors, one per entry of `cols`, each of the
recycled length, named by `names(cols)`.

## The closed form the kernel evaluates

With \\c = 1 + 2\sigma\mu\\ and \\\alpha = \sqrt{c}/\sigma\\, the
half-integer order collapses the Bessel function to a finite sum and the
log-likelihood to \$\$\ell(y) = y\log\mu - \tfrac{y}{2}\log c +
\tfrac{1}{\sigma} + \psi(\alpha) - \log y!,\qquad \psi(\alpha) =
-\alpha + \log S_y(\alpha),\$\$ \$\$S_y(\alpha) = \sum\_{k=0}^{y-1}
\dfrac{\Gamma(y+k)}{\Gamma(k+1)\Gamma(y-k)}\\(2\alpha)^{-k}, \qquad S_0
= 1.\$\$ \\S_y\\ sums \\y\\ positive terms on the log scale, so nothing
cancels. The derivatives of \\\psi\\ in \\\alpha\\ are the weighted
rising-factorial moments of \\k\\ under those terms; everything else is
elementary.

## Which kernel runs

`pig1_hd_cpp` and `pig2_hd_cpp` are **explicit closed-form kernels**,
every partial written out by hand, and they are what the package methods
run. A second pair, `pig1_hd_jet_cpp` and `pig2_hd_jet_cpp`, carries a
bivariate jet truncated at total order four through the same expression;
it shares no algebra with the explicit route, so the tests compare the
two with no tolerance to hide behind, and it is not used in production.
Measured over \\2\times10^5\\ observations, the explicit kernel takes
0.24 seconds against the jet's 1.30, and the two agree to
\\5\times10^{-14}\\.

## Rows that are not admissible

A `y` that is negative, non-integer or non-finite never reaches the
kernel and gets a row of `NaN`.
[`distrib_pdf.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Pig1Distrib.md)
turns that into a log-probability of `-Inf`, the right answer for a
point off the support.

## Notation

\\\mu\\ is the mean, \\\sigma\\ the dispersion, \\\alpha\\ the Bessel
argument, \\K\_\nu\\ the modified Bessel function of the second kind,
and \\\ell\\ the log-mass of one observation.

## See also

[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
and
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
for the two families, and
[`distrib_gradient.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Pig1Distrib.md)
for a method that calls this.

## Examples

``` r
# The score, taken directly.
distributions7:::pig_hd_block(0:4, list(mu = 3, sigma = 0.8),
                              c(mu = "d10", sigma = "d01"),
                              distributions7:::pig1_hd_cpp)
#> $mu
#> [1] -0.41522740 -0.21982510 -0.05881615  0.07269097  0.18188693
#> 
#> $sigma
#> [1]  0.6433956  0.1261542 -0.2083725 -0.3861708 -0.4454410
#> 

# ...which is what the method returns.
all.equal(distributions7:::pig_hd_block(0:4, list(mu = 3, sigma = 0.8),
                                        c(mu = "d10", sigma = "d01"),
                                        distributions7:::pig1_hd_cpp),
          distrib_gradient(pig1_distrib(), 0:4,
                           list(mu = 3, sigma = 0.8)))
#> [1] TRUE

# The explicit kernel and the jet twin agree, sharing no algebra.
y <- as.numeric(0:6)
max(abs(distributions7:::pig1_hd_cpp(y, rep(3, 7), rep(0.8, 7), 1L) -
        distributions7:::pig1_hd_jet_cpp(y, rep(3, 7), rep(0.8, 7))))
#> [1] 5.151435e-14

# A point off the support gives NaN here and -Inf in the density.
distributions7:::pig_hd_block(c(-1, 1.5, 2), list(mu = 3, sigma = 0.8),
                              c(l = "l"), distributions7:::pig1_hd_cpp)
#> $l
#> [1]       NaN       NaN -1.727361
#> 
```
