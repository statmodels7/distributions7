# Standard Deviations and Correlations of a Multivariate Gaussian

The standard deviations and correlations of the response, whichever side
the structure parametrises. A precision structure additionally reports
the **partial** correlations, which are what it describes directly:
\\-\Omega\_{jk}/\sqrt{\Omega\_{jj}\Omega\_{kk}}\\ is the correlation of
two coordinates given all the others, and it is zero exactly where the
precision has a zero.

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- theta:

  A named list of parameters.

- ...:

  Unused.

## Value

A list as described in
[`mv_derived`](https://statmodels7.github.io/distributions7/reference/mv_derived.md).
