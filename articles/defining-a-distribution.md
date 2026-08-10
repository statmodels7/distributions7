# Defining a new distribution

``` r

library(distributions7)
```

## The idea in one sentence

A new distribution in **distributions7** requires only its **density**
(`distrib_pdf`). Everything else — cumulative distribution function,
quantile function, random number generator, score (gradient), observed
and expected Hessian, and the moments — is then available automatically
through numerical fallbacks. A closed-form expression derived later is
registered as a method and transparently takes over, with no change to
calling code.

## How dispatch makes this work

Every operation in the package is an S7 *generic* dispatching on the
distribution object:

| Generic | Meaning |
|----|----|
| `distrib_pdf` | density / mass function |
| `distrib_cdf` | cumulative distribution function |
| `distrib_quantile` | quantile function (inverse CDF) |
| `distrib_rng` | random number generator |
| `distrib_gradient` | score: first derivatives of the log-density w.r.t. the parameters |
| `distrib_hessian` | observed Hessian of the log-density |
| `distrib_expected_hessian` | expected Hessian (negative Fisher information) |
| `distrib_deriv3`, `distrib_deriv4` | third- and fourth-order derivatives (observed or expected) |
| `distrib_grad_y`, `distrib_hess_y` | derivatives of the log-density w.r.t. the response `y` (continuous) |
| `mean`, `variance`, `std_dev`, `skewness`, `kurtosis`, `moment` | moments |

Default methods for all of these (except `distrib_pdf`) are registered
on the **base classes** `continuous_distrib` and `discrete_distrib`.
Because S7 always selects the most specific applicable method, a
subclass inherits those defaults until something more specific is
registered. The defaults are:

- `distrib_cdf` — numerical integration (continuous) or summation of the
  pmf (discrete) of the density;
- `distrib_quantile` — root-finding / table inversion on the CDF;
- `distrib_rng` — Generalized Ratio-of-Uniforms
  ([`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)),
  which needs nothing but the density; inverse-transform sampling,
  `distrib_quantile(runif(n))`, is used instead whenever an analytical
  quantile function is supplied;
- `distrib_gradient`, `distrib_hessian` — central finite differences of
  the log-density (see
  [`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md)
  /
  [`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md));
- `distrib_deriv3`, `distrib_deriv4` — finite differences of the Hessian
  (see
  [`numerical_deriv3()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv3.md)
  /
  [`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)),
  with the expected versions obtained by
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
  of the observed ones;
- `distrib_expected_hessian` — the
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
  of the observed Hessian;
- the moments — numerical
  [`expectation()`](https://statmodels7.github.io/distributions7/reference/expectation.md)
  of the relevant power.

## The minimal recipe

A new distribution needs exactly three things:

1.  **A class** that inherits from `continuous_distrib` or
    `discrete_distrib`.
2.  **A `distrib_pdf` method** for that class.
3.  **A constructor** that fills in the metadata (name, support bounds,
    parameter names, domains and link functions).

Nothing else is required.

### Worked example: the Laplace distribution

The Laplace (double-exponential) density is
$`f(y;\mu,b)=\frac{1}{2b}\exp\!\left(-\frac{|y-\mu|}{b}\right)`$, with
$`\mu\in\mathbb{R}`$ and $`b>0`$. Base R has no `dlaplace`, so this is a
genuinely new distribution.

**Step 1 — the class.** It is continuous, so we inherit from
`continuous_distrib`:

``` r

Laplace <- S7::new_class("Laplace", parent = continuous_distrib)
```

**Step 2 — the density.** A method for `distrib_pdf`. The signature is
always `(distrib, y, theta, log = FALSE)`, where `theta` is a named list
of parameters (each element a scalar or a vector aligned with `y`):

``` r

S7::method(distrib_pdf, Laplace) <- function(distrib, y, theta, log = FALSE) {
  mu <- theta[[1]]
  b  <- theta[[2]]
  log_d <- -log(2 * b) - abs(y - mu) / b
  if (log) log_d else exp(log_d)
}
```

Two conventions matter. `theta` is indexed positionally (`theta[[1]]`,
`theta[[2]]`) so that the parameter *names* stay free, and the `log`
branch is always implemented, because the numerical machinery evaluates
the density on the log scale for stability.

**Step 3 — the constructor.** A plain function returning an object of
the class with the metadata filled in:

``` r

laplace_distrib <- function(link_mu = linkfunctions7::identity_link(),
                            link_b  = linkfunctions7::log_link()) {
  Laplace(
    distrib_name = "laplace",
    dimension    = "univariate",
    bounds       = c(-Inf, Inf),                 # support of Y
    params       = c("mu", "b"),
    params_interpretation = c(mu = "location", b = "scale"),
    n_params     = 2,
    params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
    link_params   = list(mu = link_mu, b = link_b)
  )
}
```

The `params_bounds` are used by the fallbacks (finite-difference steps
are shrunk so they never step outside a parameter’s domain) and by
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md).
The `link_params` are used by downstream modeling code;
`identity_link()` and `log_link()` come from **linkfunctions7**.

That is the entire definition. Everything now works:

``` r

d  <- laplace_distrib()
th <- list(mu = 1, b = 2)

# density (what we defined)
distrib_pdf(d, c(0, 1, 3), th)
#> [1] 0.15163266 0.25000000 0.09196986

# CDF, quantile, RNG --- all numerical fallbacks
distrib_cdf(d, c(-1, 1, 4), th)
#> [1] 0.1839397 0.5000000 0.8884349
distrib_quantile(d, c(0.25, 0.5, 0.75), th)   # median is exactly mu = 1
#> [1] -0.3862944  1.0000000  2.3862944
head(distrib_rng(d, 5, th))
#> [1]  5.723932 -2.857153  0.855691  2.715865  0.144959

# score and Hessian --- finite differences of the log-density
str(distrib_gradient(d, c(-1, 3), th))
#> List of 2
#>  $ mu: num [1:2] -0.5 0.5
#>  $ b : num [1:2] 1.83e-11 1.83e-11
str(distrib_hessian(d, c(-1, 3), th))
#> List of 3
#>  $ mu_mu: num [1:2] 0 0
#>  $ b_b  : num [1:2] -0.25 -0.25
#>  $ mu_b : num [1:2] 0.25 -0.25

# moments --- numerical integration against the density
c(mean = mean(d, th), variance = variance(d, th))  # Laplace variance is 2 b^2 = 8
#>     mean variance 
#>        1        8
```

We can check the fallbacks against the known closed forms for the
Laplace:

``` r

# analytical Laplace CDF
lap_cdf <- function(q, mu, b) {
  ifelse(q < mu, 0.5 * exp((q - mu) / b), 1 - 0.5 * exp(-(q - mu) / b))
}
q <- c(-3, -1, 0, 1, 4)
max(abs(distrib_cdf(d, q, th) - lap_cdf(q, 1, 2)))   # ~1e-10
#> [1] 1.387779e-15
```

``` r

# variance should be exactly 2 b^2
all.equal(variance(d, th), 2 * th$b^2)
#> [1] TRUE
```

The object also prints and plots like any built-in distribution:

``` r

d
#> Distribution: Laplace
#> Type:         Continuous
#> Dimensions:   univariate
#> 
#> Parameters:
#>   mu (location)           | Link: identity   | Domain: (-Inf, Inf)
#>   b  (scale)              | Link: log        | Domain: (0, Inf)
```

``` r

plot(d, th)
```

![](defining-a-distribution_files/figure-html/unnamed-chunk-9-1.png)

## Discrete distributions

For a discrete distribution, inherit from `discrete_distrib` instead and
let `distrib_pdf` be the **probability mass function**. The one extra
requirement is that the support has a **finite lower bound**
(`bounds[1]`), which every standard count distribution satisfies. The
discrete fallbacks build a cumulative-mass table from that lower bound,
growing it geometrically as needed.

``` r

# A shifted-geometric example: P(Y = k) = (1-p) p^k for k = 0, 1, 2, ...
Geom <- S7::new_class("Geom", parent = discrete_distrib)
S7::method(distrib_pdf, Geom) <- function(distrib, y, theta, log = FALSE) {
  p <- theta[[1]]
  log_d <- log(1 - p) + y * log(p)
  if (log) log_d else exp(log_d)
}
geom_distrib <- function(link_p = linkfunctions7::logit_link()) {
  Geom(
    distrib_name = "geometric", dimension = "univariate", bounds = c(0, Inf),
    params = "p", params_interpretation = c(p = "success prob."),
    n_params = 1, params_bounds = list(p = c(0, 1)),
    link_params = list(p = link_p)
  )
}

g <- geom_distrib()
distrib_cdf(g, 0:3, list(p = 0.6))                      # 1 - p^(k+1)
#> [1] 0.4000 0.6400 0.7840 0.8704
distrib_quantile(g, c(0.5, 0.9), list(p = 0.6))
#> [1] 1 4
c(mean = mean(g, list(p = 0.6)))                        # p / (1 - p) = 1.5
#> mean 
#>  1.5
```

## Adding closed forms for performance

The numerical fallbacks are convenient but slower and less precise than
analytical formulas. Whenever a closed form is available, registering it
as a method makes it take over automatically, and existing code keeps
working. There is no need to override everything: any subset of the
methods can be supplied.

For the Laplace, all of these are available in closed form. For instance
the CDF, quantile and score:

``` r

S7::method(distrib_cdf, Laplace) <- function(distrib, q, theta,
                                             lower.tail = TRUE, log.p = FALSE) {
  mu <- theta[[1]]; b <- theta[[2]]
  res <- ifelse(q < mu, 0.5 * exp((q - mu) / b), 1 - 0.5 * exp(-(q - mu) / b))
  if (!lower.tail) res <- 1 - res
  if (log.p) log(res) else res
}

S7::method(distrib_quantile, Laplace) <- function(distrib, p, theta,
                                                  lower.tail = TRUE, log.p = FALSE) {
  mu <- theta[[1]]; b <- theta[[2]]
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  mu - b * sign(p - 0.5) * log(1 - 2 * abs(p - 0.5))
}

# score: d/dmu = sign(y - mu) / b,  d/db = (|y - mu| / b - 1) / b
S7::method(distrib_gradient, Laplace) <- function(distrib, y, theta,
                                                  scale = c("parameter", "link"), ...) {
  mu <- theta[[1]]; b <- theta[[2]]
  list(
    mu = sign(y - mu) / b,
    b  = (abs(y - mu) / b - 1) / b
  )
}
```

Note the signature of the derivative method: it must mirror the generic,
which carries a `scale` argument (and `...`). **A method always returns
parameter-scale derivatives**; the chain rule is never implemented by
hand. The generic applies the link-scale transformation afterwards when
the caller asks for `scale = "link"`. The same applies to
`distrib_hessian`, `distrib_expected_hessian`, `distrib_deriv3` and
`distrib_deriv4`; the probability-function methods (`distrib_pdf`,
`distrib_cdf`, `distrib_quantile`, `distrib_rng`) keep their plain
signatures.

The results are identical to what the fallbacks produced, only faster
and exact:

``` r

max(abs(distrib_cdf(d, q, th) - lap_cdf(q, 1, 2)))       # now exact
#> [1] 0
str(distrib_gradient(d, c(-1, 3), th))                   # now analytical
#> List of 2
#>  $ mu: num [1:2] -0.5 0.5
#>  $ b : num [1:2] 0 0
```

A practical suggestion for the order in which to add closed forms:

1.  **`distrib_gradient`** — the single biggest speed-up for model
    fitting, since the score is evaluated at every optimization step.
    Writing it in C++ (via **Rcpp**), as the built-in distributions do,
    helps further.
2.  **`distrib_hessian`** — for Newton-type optimization and standard
    errors.
3.  **`distrib_cdf` / `distrib_quantile` / `distrib_rng`** — for fast or
    high-volume simulation and tail probabilities.
4.  **`distrib_expected_hessian`** and **analytical moments** — when a
    closed form exists and the fitting method is Fisher scoring, or
    exact moments are needed.

## Checking analytical derivatives

Hand-written analytical derivatives should be verified against the
numerical reference implementations, which are exported for exactly this
purpose:

``` r

y <- distrib_rng(d, 50, th)

ana <- distrib_gradient(d, y, th)
num <- numerical_gradient(d, y, th)
max(abs(ana$mu - num$mu))
#> [1] 2.489187e-11
max(abs(ana$b  - num$b))
#> [1] 5.887246e-11
```

[`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md)
and
[`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md)
take the same `(distrib, y, theta)` arguments as the generics and return
the finite-difference estimates, so a component-by-component comparison
is a one-liner.

## Validating the distribution

Once a distribution is defined,
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
runs a battery of numerical self-consistency checks: it verifies that
the density integrates (or sums) to one, that the CDF, quantile function
and RNG agree with each other, and that every analytical derivative
matches finite differences.

``` r

check_distrib(laplace_distrib(), theta = list(mu = 1, b = 2),
              n = 40, nsim = 2e4, orders = 1:2)
#> Distribution: laplace
#> Parameters:   mu = 1, b = 2
#> Observations: 40   Monte Carlo: 20000
#> 
#>   [OK  ] density integrates to 1                     1.81e-10
#>   [OK  ] density is non-negative                     5.00e-03
#>   [OK  ] cdf in [0,1] and non-decreasing             2.46e-02
#>   [OK  ] cdf agrees with the density                 6.25e-07
#>   [OK  ] quantile/cdf round-trip                     1.39e-17
#>   [OK  ] rng matches the cdf                         1.70e+00
#>   [OK  ] gradient vs finite differences              6.37e-11
#>   [OK  ] hessian vs finite differences               0.00e+00
#>   [OK  ] expected information vs Monte Carlo         9.33e-01
#>   [OK  ] response derivatives vs finite differences  0.00e+00
#>   [OK  ] link-scale gradient vs finite differences   3.20e-09
#> 
#> All 11 checks passed.
```

This is the fastest way to catch a mistake in a hand-derived score or
Hessian: an analytical `distrib_gradient` that disagrees with the
density makes the corresponding row come back `FAIL`.

## Derivatives on the link scale

Every parameter carries a link function, and the derivative generics
accept `scale = "link"` to differentiate with respect to the
unconstrained linear predictor $`\eta = g(\theta)`$ instead of
$`\theta`$ itself:

``` r

d <- gaussian1_distrib()          # identity link on mu, log link on sigma
th <- list(mu = 1.5, sigma = 2)
y <- c(0.4, 1.9, 3.2)

distrib_gradient(d, y, th)                    # d l / d theta
#> $mu
#> [1] -0.275  0.100  0.425
#> 
#> $sigma
#> [1] -0.34875 -0.48000 -0.13875
distrib_gradient(d, y, th, scale = "link")    # d l / d eta
#> $mu
#> [1] -0.275  0.100  0.425
#> 
#> $sigma
#> [1] -0.6975 -0.9600 -0.2775
```

Because $`\sigma = e^{\eta}`$, the second element of the link-scale
score is simply the parameter-scale score times $`\sigma`$. Higher
orders follow the Faà di Bruno chain rule and are available up to order
4; see
[`?link_scale_derivatives`](https://statmodels7.github.io/distributions7/reference/link_scale_derivatives.md)
for the formulas.

## Fitting a distribution to data

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
maximizes the likelihood on the link scale — where the parameters are
unconstrained — using the analytical score and information, and then
reports the estimates back on the parameter scale:

``` r

set.seed(1)
d <- gaussian1_distrib()
y <- distrib_rng(d, 500, list(mu = 2, sigma = 3))
fit_distrib(d, y)
#> Maximum-likelihood fit: gaussian1
#> Observations: 500   Log-likelihood: -1264   AIC: 2532   BIC: 2541
#> Method: Fisher scoring   iterations: 5   evaluations: f 6, g 6   time: 17 ms
#> Converged: yes (gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative))
#> 
#> Parameter scale:
#>       Estimate Std. Error   2.5%  97.5%
#> mu      2.0679     0.1356 1.8021 2.3338
#> sigma   3.0327     0.0959 2.8505 3.2267
#> 
#> Link scale:
#>       Estimate Std. Error   2.5%  97.5%
#> mu      2.0679     0.1356 1.8021 2.3338
#> sigma   1.1095     0.0316 1.0475 1.1714
```

Confidence intervals are built symmetrically on the link scale and
mapped through $`g^{-1}`$, so they always respect the parameter’s
domain. For a probability, for instance, the interval can never leave
$`(0, 1)`$:

``` r

set.seed(2)
fit_distrib(bernoulli_distrib(), rbinom(50, 1, 0.9))
#> Maximum-likelihood fit: bernoulli
#> Observations: 50   Log-likelihood: -20.25   AIC: 42.5   BIC: 44.41
#> Method: Fisher scoring   iterations: 6   evaluations: f 10, g 7   time: 4 ms
#> Converged: yes (gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative))
#> 
#> Parameter scale:
#>    Estimate Std. Error   2.5%  97.5%
#> mu     0.86     0.0491 0.7343 0.9318
#> 
#> Link scale:
#>    Estimate Std. Error   2.5%  97.5%
#> mu   1.8153     0.4076 1.0165 2.6141
```

## Non-differentiable parameters

Some distributions have a log-likelihood that is not differentiable in
one of its parameters. The archetype is the **Laplace** distribution
$`f(y;\mu,b)=\frac{1}{2b}e^{-|y-\mu|/b}`$, whose density has a kink at
$`y=\mu`$: the score
$`\partial\ell/\partial\mu=\operatorname{sign}(y-\mu)/b`$ exists only
almost everywhere, and the observed curvature
$`\partial^2\ell/\partial\mu^2`$ is **zero**, so Newton-Raphson cannot
update $`\mu`$.

The package handles this in two ways. First, a distribution can declare
which parameters are smooth via the `params_smooth` metadata (shown in
[`print()`](https://rdrr.io/r/base/print.html)):

``` r

d <- laplace_distrib()
param_smoothness(d)
#>   mu    b 
#> TRUE TRUE
```

Second, and most importantly, `distrib_expected_hessian` returns the
**Fisher information** computed from the variance of the score,
$`I=\mathbb{E}[\nabla\ell\,\nabla\ell^\top]`$, rather than from
$`-\mathbb{E}[H]`$. For a regular model the two coincide, by the second
Bartlett identity, but only the score-variance form stays valid when the
identity fails. It gives the correct information $`1/b^2`$ for the
Laplace location, where $`-\mathbb{E}[H]`$ would give a degenerate zero:

``` r

th <- list(mu = 1, b = 2)
distrib_hessian(d, c(-1, 0, 2), th)$mu_mu          # observed curvature: 0
#> [1] 0 0 0
distrib_expected_hessian(d, 0, th)$mu_mu           # Fisher information: -1/b^2
#> [1] -0.25
```

This means that even a user-defined non-smooth distribution supplying
only its density (and, ideally, an analytic score at the kink) obtains
the right expected information automatically, because the default
`distrib_expected_hessian` uses the outer product of the score. See
[`?laplace_distrib`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
for the full worked example.

## Summary

- **Minimum to define a distribution:** a subclass of
  `continuous_distrib` or `discrete_distrib`, a `distrib_pdf` method,
  and a constructor with the metadata. (Discrete distributions
  additionally need a finite lower support bound.)
- The density alone yields the CDF, quantile function, RNG, score,
  observed and expected Hessian, and all moments, via numerical
  fallbacks.
- Register analytical methods incrementally for speed and precision;
  they override the fallbacks automatically through S7 dispatch.
- Validate analytical derivatives with
  [`numerical_gradient()`](https://statmodels7.github.io/distributions7/reference/numerical_gradient.md)
  /
  [`numerical_hessian()`](https://statmodels7.github.io/distributions7/reference/numerical_hessian.md).
