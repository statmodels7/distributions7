# Register All Four CDF Orders on a Sum of Normal Tails

Turns one term function into the four S7 methods a family of the shape
\\F = c_0 + \sum_k s_k e^{w_k}\Phi(x_k)\\ needs, so that the family
states its terms once instead of four times. The inverse Gaussian and
the elastic net are the two families registered through it.

## Usage

``` r
register_phi_terms_cdf(cls, term_fn)
```

## Arguments

- cls:

  The S7 class to register on.

- term_fn:

  A function of `(distrib, q, theta)` returning the list of terms
  [`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md)
  documents.

## Value

Invisibly `NULL`. Called for the registration.

## Details

All four orders are registered,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
included, so these families take the closed route from the first order
up and never reach a stencil. `force(o)` inside the factory is what
keeps the four registrations from sharing one order.

## See also

[`phi_terms_cdf_deriv_k()`](https://statmodels7.github.io/distributions7/reference/phi_terms_cdf_deriv_k.md),
the body it registers;
[`distrib_grad_cdf.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.InvGauss1Distrib.md)
and
[`distrib_grad_cdf.EnetDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.EnetDistrib.md),
the two families.
