# The Pseudo-Huber Does Not Write Its Expected Information Out

The method above calls
[`expected_derivative`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
and then replaces two components that vanish by symmetry, so the
registration improves the approximation rather than replacing it.

## Arguments

- x:

  A `PseudoHuberDistrib` object.

- ...:

  Unused.

## Value

`FALSE`.

## Details

Measured at 100 observations it costs 10980 ms, where the families that
write their expected information out answer in a median of 0.183 ms.
Reported as exact, it made
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reject a legitimate `fisher_scoring(approx = )` here with a message that
was untrue, and take a quadrature it believed was a formula when
assembling the standard errors.

## See also

[`expected_hessian_exact`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
