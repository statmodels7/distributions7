# Student t Distribution in the Standard Deviation, Obtained

The same family as
[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md),
obtained through
[`reparametrize()`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
on
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
rather than written out.

## Usage

``` r
student_t2_by_reparam(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_nu = bounded_link(lwr = 2)
)
```

## Arguments

- link_mu:

  Link function for the location. Defaults to the identity.

- link_sigma:

  Link function for the standard deviation. Defaults to the log.

- link_nu:

  Link function for the degrees of freedom. Defaults to a link bounded
  below at two.

## Value

A reparametrized distribution object.

## Details

This exists as a reference for the tests and is not exported:
[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md)
carries its own kernels, so the two are independent implementations of
one law.

## See also

[`student_t2_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t2_distrib.md)
