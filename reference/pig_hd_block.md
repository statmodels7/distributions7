# Log-Likelihood Derivatives of the Poisson-Inverse Gaussian

Evaluates the log-likelihood and its fourteen partial derivatives to
fourth order at once, through the compiled kernel, and returns the block
of components an order asks for.

## Usage

``` r
pig_hd_block(y, theta, cols, kernel, threads = 1L)
```

## Arguments

- y:

  A numeric vector of observations.

- theta:

  The aligned parameter list.

- cols:

  The kernel columns wanted, by name.

- kernel:

  The compiled kernel, `pig1_hd_cpp` or `pig2_hd_cpp`.

- threads:

  How many threads that kernel may use.

## Value

A named list of derivative component vectors.

## Details

With \\c = 1 + 2\sigma\mu\\ and \\\alpha = \sqrt{c}/\sigma\\, the
half-integer order collapses the Bessel function to a finite sum and the
log-likelihood to \$\$\ell(y) = y\log\mu - \tfrac{y}{2}\log c +
\tfrac{1}{\sigma} + \psi(\alpha) - \log y!,\qquad \psi(\alpha) =
-\alpha + \log S_y(\alpha),\$\$ where \\S_y\\ sums \\y\\ positive terms
on the log scale. The kernel carries a bivariate jet truncated at total
order four through this expression, so every partial is exact and no
chain rule is transcribed.
