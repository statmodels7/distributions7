# Validate a Multivariate Distribution

Runs the battery of
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
on a multivariate distribution, where the one-dimensional checks have no
counterpart.

## Usage

``` r
check_distrib_mv(distrib, theta, n, nsim, tol)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters, already aligned.

- n:

  The number of observations drawn for the derivative checks.

- nsim:

  The Monte Carlo sample size.

- tol:

  The tolerance of the gradient and response-derivative comparisons,
  applied to the absolute difference divided by the larger of one and
  the magnitude of the reference; the Hessian comparison uses
  `10 * tol`. The other checks use fixed thresholds.

## Value

A list of one-row data frames, as
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
builds.

## Details

Five of the univariate checks have no counterpart in \\p\\ dimensions:
the integration of the density by one-dimensional quadrature, the two
checks of the distribution function (its range and its agreement with
the density), the quantile round trip, and the comparison of the
generator with the distribution function. The distribution function is
an integral over an orthant, and the quantile function inverts an
ordering that does not exist in several dimensions. Reporting these
checks as failures would mark a correct family as defective, and a user
validating their own distribution could not tell a real defect from
them.

The battery runs the checks that do generalize:

1.  **the density is positive and finite** on a sample;

2.  **the density integrates to one**. For a family that enumerates its
    support through
    [`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
    this is an exact sum over that support, with tolerance `1e-10`;
    otherwise it is importance sampling from the proposal that
    [`mv_reference_draw()`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
    supplies, with tolerance `20 / sqrt(nsim)`. The default proposal is
    a gaussian with the same mean and twice the covariance; the Student
    \\t\\ uses a Student \\t\\ with half the degrees of freedom and
    twice the scale matrix, and the Dirichlet the uniform law on the
    simplex. The proposal is not the distribution itself, which would
    make the ratio identically one;

3.  **the score has mean zero**, the first Bartlett identity, under the
    distribution's own generator, so it is also a check that the
    generator and the density describe the same law;

4.  **gradient and Hessian** against finite differences of the summed
    log-density;

5.  **the expected information** against the Monte Carlo average of the
    observed Hessian, and against the variance of the score, which is
    the second Bartlett identity;

6.  **the generator** against the first two moments. When the mean or
    the covariance of the family is not finite (a Student \\t\\ with
    \\\nu \le 2\\) the row has no statistic, and
    [`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
    lists it among the skipped checks with that reason. The tolerance is
    `20 / sqrt(nsim)`, which assumes a finite fourth moment: for a
    Student \\t\\ with \\2 \< \nu \le 4\\ the sample covariance has no
    variance, and the row can fail with a correct generator;

7.  **the response derivatives** against finite differences in \\y\\.

The last of these is run only when it applies, as the univariate battery
omits the checks that have no counterpart for a discrete family. A
family with an enumerable support is discrete and has no derivative in
the response, and the multivariate base class signals an error for
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md),
so a family that has not registered a method has no response derivatives
by design, and the check is omitted for it.

## See also

[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
