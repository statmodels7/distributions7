# The Centered Skew Normal Does Not Write Its Expected Information Out

The method above is the CHAIN onto
[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
whose expected information is the base class's quadrature, so the
registration says where the arithmetic is assembled and not that it is
closed form.

## Arguments

- x:

  A `SkewNormal2Distrib` object.

- ...:

  Unused.

## Value

`FALSE`.

## Details

Measured at 100 observations it costs 5220 ms against the parent's 2230
– more than what it chains onto, the chain being paid on top – where the
families that do write it out answer in a median of 0.183 ms. Reported
as exact, it made
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reject a legitimate `fisher_scoring(approx = )` here with a message that
was untrue.

## See also

[`expected_hessian_exact`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.md)
