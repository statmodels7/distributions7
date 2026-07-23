# Calculate the Expected Value of a Function

Computes the expected value of a given function \\f(y)\\ with respect to
a probability distribution defined by `distrib`. It automatically
handles continuous distributions (via numerical integration) and
discrete distributions (via series summation).

## Usage

``` r
expectation(distrib, f, theta, ...)
```

## Arguments

- distrib:

  An object of class `"distrib"`

- f:

  A function representing the transformation of the random variable
  \\y\\. **Signature:** It must accept arguments `y`, `theta`, and `...`
  (see Details).

- theta:

  A named list of parameters for the distribution (e.g.,
  `list(mu=10, sigma=2)`). Vectors inside this list allow computing
  expectations for multiple distribution parametrizations at once.

- ...:

  Additional arguments passed directly to the function `f`.
  **Vectorization:** These arguments are fully vectorized. If vectors
  are provided, they are recycled against the parameters in `theta`
  according to standard R recycling rules.

## Value

A numeric vector containing the expected values. The length corresponds
to the maximum length among all vectors in `theta` and `...`.

## Details

The function calculates:

- \\E\[f(Y)\] = \int\_{lb}^{ub} f(y, \theta, \dots) \cdot p(y\|\theta)
  \\ dy\\ (Continuous)

- \\E\[f(Y)\] = \sum\_{y=lb}^{ub} f(y, \theta, \dots) \cdot
  P(y\|\theta)\\ (Discrete)

**Vectorization:** The function iterates over the longest vector found
among `theta` and `...`. For example, if `theta$mu` has length 2 and you
pass a vector of length 2 to `...`, the function computes the
expectation for the paired values. If lengths differ, standard R
recycling applies.

**Requirements for `f`:** The user-provided function `f` must be defined
with the signature: `f(y, theta, ...)`

## Examples

``` r
if (FALSE) { # \dontrun{
distrib <- poisson_distrib()

# Define f accepting y, theta, and extra parameter gamma
f_pow <- function(y, theta, gamma = 1) {
  y^gamma
}

# --- Example 1: Basic usage ---
expectation(distrib, f_pow, theta = list(mu = 10), gamma = 2)
} # }
```
