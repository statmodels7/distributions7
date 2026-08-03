# Derivatives and the link scale

``` r

library(distributions7)
```

Most of what **distributions7** offers beyond densities is derivatives:
the score, the information, and the third and fourth derivatives of the
log-likelihood. This vignette explains what they are, the two scales
they are available on, and why the second scale — the link scale — is
the one an optimiser needs.

## The derivatives of the log-likelihood

For one observation, write $`\ell(\theta) = \log f(y; \theta)`$. The
package returns its derivatives with respect to the parameters, per
observation, as named lists:

``` r

d <- gaussian_distrib()
theta <- list(mu = 2, sigma = 3)
y <- distrib_rng(d, 4, theta)

distrib_gradient(d, y, theta)          # score: d ell / d theta
#> $mu
#> [1] -0.20881794  0.06121444 -0.27854287  0.53176027
#> 
#> $sigma
#> [1] -0.2025185 -0.3220917 -0.1005749  0.5149736
```

``` r

distrib_hessian(d, y, theta)           # observed: d^2 ell / d theta^2
#> $mu_mu
#> [1] -0.1111111 -0.1111111 -0.1111111 -0.1111111
#> 
#> $sigma_sigma
#> [1] -0.01970368  0.09986949 -0.12164728 -0.73719583
#> 
#> $mu_sigma
#> [1]  0.13921196 -0.04080963  0.18569525 -0.35450684
```

The **expected** Hessian — the negative Fisher information — is a
separate generic, because for a correct model it does not depend on
$`y`$:

``` r

distrib_expected_hessian(d, 0, theta)
#> $mu_mu
#> [1] -0.1111111
#> 
#> $sigma_sigma
#> [1] -0.2222222
#> 
#> $mu_sigma
#> [1] 0
```

Third and fourth derivatives are there too, keyed by their multi-index,
for anything that needs curvature beyond the quadratic:

``` r

distrib_deriv3(d, y, theta)
#> $mu_mu_mu
#> [1] 0 0 0 0
#> 
#> $mu_mu_sigma
#> [1] 0.07407407 0.07407407 0.07407407 0.07407407
#> 
#> $mu_sigma_sigma
#> [1] -0.13921196  0.04080963 -0.18569525  0.35450684
#> 
#> $sigma_sigma_sigma
#> [1]  0.10034565 -0.05908524  0.23627045  1.05700185
```

These are exact, computed by closed-form C++ kernels, not by
differencing. A separate pair,
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
and
[`distrib_hess_y()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.md),
differentiates with respect to the response $`y`$ instead of the
parameters.

## Two scales

Every derivative generic takes a `scale` argument. The default,
`scale = "parameter"`, is what the examples above used: derivatives with
respect to $`\theta`$ on its natural, constrained scale.

`scale = "link"` instead gives derivatives with respect to the
**unconstrained** parameter $`\eta = g(\theta)`$ behind each link. For
the Gaussian the default link on $`\sigma`$ is the log, so
$`\eta_\sigma = \log\sigma`$:

``` r

distrib_gradient(d, y, theta, scale = "link")
#> $mu
#> [1] -0.20881794  0.06121444 -0.27854287  0.53176027
#> 
#> $sigma
#> [1] -0.6075556 -0.9662751 -0.3017248  1.5449208
```

Nothing was differentiated numerically to get this. The chain rule was
applied to the parameter-scale derivatives, exactly. For a single
parameter with link $`g`$,
$`\partial\ell/\partial\eta = (\partial\ell/\partial\theta)\,/\,g'(\theta)`$;
with the log link on $`\sigma`$ that factor is simply $`\sigma`$:

``` r

g_par  <- distrib_gradient(d, y, theta)
g_link <- distrib_gradient(d, y, theta, scale = "link")

all.equal(g_par$sigma * theta$sigma, g_link$sigma)
#> [1] TRUE
```

Under an identity link the two scales coincide, as they must:

``` r

d_id <- gaussian_distrib(
  link_mu    = linkfunctions7::identity_link(),
  link_sigma = linkfunctions7::identity_link()
)
all.equal(
  distrib_gradient(d_id, y, theta),
  distrib_gradient(d_id, y, theta, scale = "link")
)
#> [1] TRUE
```

The transformation works at every order. The general rule is Faà di
Bruno’s formula — the chain rule for higher derivatives — applied with a
diagonal Jacobian, since each parameter has its own link. It is written
once, in the body of each generic, so `scale = "link"` is available
uniformly:

``` r

distrib_deriv3(d, y, theta, scale = "link")
#> $mu_mu_mu
#> [1] 0 0 0 0
#> 
#> $mu_mu_sigma
#> [1] 0.2222222 0.2222222 0.2222222 0.2222222
#> 
#> $mu_sigma_sigma
#> [1] -0.8352717  0.2448578 -1.1141715  2.1270411
#> 
#> $sigma_sigma_sigma
#> [1]  1.5697775  0.1348995  2.7931007 10.1796834
```

## Why the link scale is the useful one

A standard deviation is positive, a probability lives in $`(0, 1)`$, a
shape parameter is positive. Maximising a likelihood over such
constrained parameters means either fighting the boundary or bolting on
constraints. Reparameterise through the link and the problem becomes
unconstrained: $`\eta`$ ranges over the whole real line, and a Newton or
Fisher-scoring step needs exactly the link-scale score and information
this package computes in closed form.

That is precisely what
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
does — it optimises on the link scale and maps back — and it is why its
confidence intervals cannot leave a parameter’s domain. See
[`vignette("fitting-a-model")`](https://statmodels7.github.io/distributions7/articles/fitting-a-model.md)
for that side of the story.

The information transforms consistently. Write $`h = g^{-1}`$, so that
$`h'(\eta) = 1/g'(\theta)`$ is the derivative of the inverse link — the
same factor that already appeared in the first-order rule above. Because
the score has zero expectation, the expected information on the link
scale is the sandwich
$`\operatorname{diag}(h') \, \mathcal{I} \, \operatorname{diag}(h')`$;
for the log link on $`\sigma`$ we have $`h'(\eta) = e^{\eta} = \sigma`$,
so the $`\sigma\sigma`$ entry is multiplied by $`\sigma^2`$:

``` r

eh_par  <- distrib_expected_hessian(d, 0, theta)
eh_link <- distrib_expected_hessian(d, 0, theta, scale = "link")

c(parameter = eh_par$sigma_sigma[1],
  scaled    = theta$sigma^2 * eh_par$sigma_sigma[1],
  link      = eh_link$sigma_sigma[1])
#>  parameter     scaled       link 
#> -0.2222222 -2.0000000 -2.0000000
```

## Expected higher derivatives, when there is no formula

[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
also take `expected = TRUE`. Some distributions have closed-form
expected derivatives; where they do not, the value is approximated, and
`approx` chooses how:

``` r

distrib_deriv3(gamma_distrib(), 0, list(mu = 3, sigma2 = 2), expected = TRUE)
#> $mu_mu_mu
#> [1] -0.2431584
#> 
#> $mu_mu_sigma2
#> [1] 0.3016318
#> 
#> $mu_sigma2_sigma2
#> [1] -0.1728947
#> 
#> $sigma2_sigma2_sigma2
#> [1] 0.2638418
```

The strategies are `"bartlett"`, `"integrate"` and `"mc"`, trading speed
against robustness. They rest on one idea: the outer product of
gradients is the order-two case of the Bartlett identities, and the
general identity expresses an expected derivative of any order through
products of lower-order ones.
[`?distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
has the details; the default is chosen to be both cheap and valid,
including for models where the usual regularity assumptions fail.

## A note for authors of derivative methods

A user-written analytical derivative method **must** include the `scale`
argument in its signature, even if the method only ever returns the
parameter scale:

``` r

S7::method(distrib_gradient, MyDist) <- function(distrib, y, theta,
                                                 scale = c("parameter", "link"), ...) {
  # ... return the parameter-scale gradient; the generic handles "link" itself
}
```

S7 requires a method’s formals to include the generic’s named arguments,
and the link-scale transformation lives in the generic, wrapping
whatever the method returns. Omit `scale` and dispatch will fail. The
numerical fallbacks and all built-in methods already follow this
pattern.

## Where to go next

- [`vignette("fitting-a-model")`](https://statmodels7.github.io/distributions7/articles/fitting-a-model.md)
  — the link scale put to work in estimation.
- [`vignette("defining-a-distribution")`](https://statmodels7.github.io/distributions7/articles/defining-a-distribution.md)
  — adding a distribution that inherits all of this.
- [`?distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
  [`?distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
  — the full reference.
