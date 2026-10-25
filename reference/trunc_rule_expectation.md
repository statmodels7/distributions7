# The Expectation of a Truncated Family on the Compiled Rule

Computes \\E_T\[f(Y)\]\\ for a
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
continuous family whose parent has a scalar route and a power-law tail,
as the sum of \\f(y_k)\\f_T(y_k)\\w_k\\ over the nodes and weights of
`trunc_cont_rule_cpp()`, the rule the compiled expected information
uses, with the anchor nodes of the tail past \\10^{60}\\ scales.

## Usage

``` r
trunc_rule_expectation(distrib, f, theta, ...)
```

## Arguments

- distrib:

  A truncated continuous family.

- f:

  The function whose expectation is taken.

- theta:

  A named list of parameters.

- ...:

  Further arguments passed to `f`.

## Value

A numeric vector with one value per parameter combination, or `NULL`.

## Details

The adaptive quadrature of
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
integrates over the whole support and reaches no tail of order
\\\|z\|^{-1-\nu}\\ at small \\\nu\\: it signals that the accuracy was
not reached, or returns a value that misses the tail. The rule
integrates the tail through the anchors, in one pass over all parameter
combinations. The rule is used for the families that
[`trunc_rule_wanted()`](https://statmodels7.github.io/distributions7/reference/trunc_rule_wanted.md)
names, the skew t and the Student t of class `StudentT1Distrib`; every
other family, a parent without a route, a constant that varies by
observation, and a rule that returns no nodes for a combination leave
the result `NULL`, and the caller takes the adaptive quadrature.

The rule is fixed, so it is accurate only for a smooth `f`. It is
therefore read by
[`trunc_expectation()`](https://statmodels7.github.io/distributions7/reference/trunc_expectation.md),
which the derivatives of the expected information call with their smooth
integrands, and not by the
[`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
method, which takes any `f` (an indicator, an absolute value) and keeps
the adaptive quadrature.
