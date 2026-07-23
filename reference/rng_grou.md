# Generalized Ratio-of-Uniforms Sampling

Draws from a continuous distribution using the Generalized
Ratio-of-Uniforms (GRoU) method. The sampler needs nothing but the (log)
density: no CDF, no quantile function and no inversion. It is the engine
behind the default
[`distrib_rng`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md)
method for continuous distributions that provide neither a native RNG
nor an analytical quantile function.

## Usage

``` r
rng_grou(distrib, n, theta, r = 2)
```

## Arguments

- distrib:

  An object inheriting from class `"continuous_distrib"`.

- n:

  Number of observations to generate.

- theta:

  A named list of parameters, each of length one.

- r:

  Numeric tuning parameter of the transformation power, default `2`.
  `r = 1` is the classical Ratio-of-Uniforms; larger values enclose
  heavier tails (`r = 2` keeps the acceptance region bounded for tails
  as heavy as the Cauchy's, which is why it is the default).

## Value

A numeric vector of length `n`.

## Details

Write \\K\\ for the density kernel. If \\(U, V)\\ is uniform over
\$\$A_r = \left\\ (u, v) : 0 \< u \leq K\\\left(v /
u^{r}\right)^{1/(r+1)} \right\\\$\$ then \\Y = V / U^{r}\\ has density
proportional to \\K\\. Sampling therefore reduces to drawing uniformly
in a bounding rectangle \\\[0, u\_{\max}\] \times \[v\_{\min},
v\_{\max}\]\\ and keeping the pairs that fall in \\A_r\\, i.e. those
with \$\$(r+1)\log U \leq \log K(Y).\$\$

Two devices make this numerically safe for an arbitrary user-supplied
density:

- **Recentring.** The kernel is shifted to its mode, \\K(z) \propto
  f(m + z)\\, and the mode is added back to the draws. Without this a
  distribution located far from the origin (say \\\mu = 1000\\) gives a
  wildly elongated bounding rectangle and an acceptance rate close to
  zero. Recentring also splits the box exactly: \\h(z) =
  z\\K(z)^{r/(r+1)}\\ is non-positive for \\z \le 0\\ and non-negative
  for \\z \ge 0\\, so \\v\_{\min}\\ and \\v\_{\max}\\ are each found on
  one side of the mode.

- **Normalisation.** The kernel is rescaled by its value at the mode, so
  that \\\max K = 1\\ and \\u\_{\max} = 1\\. Every quantity stays in a
  safe numerical range whatever the height of the density, and all
  computations are carried out on the log scale.

The bounds \\v\_{\min}, v\_{\max}\\ are obtained by expanding
geometrically away from the mode until \\h\\ has clearly turned back
towards zero (finite support bounds are used directly), then combining
[`optimize`](https://rdrr.io/r/stats/optimize.html) with a grid search
over the resulting bracket; the box is enclosing by construction for a
unimodal kernel. Candidates are generated and filtered in vectorised
batches whose size adapts to the observed acceptance rate.

## Unbounded densities

The method needs the supremum of the kernel to be attained, so a density
that diverges at a finite edge of its support would have its acceptance
region clipped. Rather than give up, such an edge is transformed away.

If \\f(y) \sim \|y-a\|^{\alpha-1}\\ near the edge \\a\\, with \\\alpha
\< 1\\, then \\X = \|Y-a\|^{1/\lambda}\\ has density proportional to
\\x^{\lambda\alpha-1}\\ there, which is bounded as soon as
\\\lambda\alpha \> 1\\. Sampling \\X\\ and mapping back is exact: no
quadrature and no inversion are involved. The exponent \\\alpha\\ is
read off the same probe that detects the divergence – walking towards
the edge in decades lifts the log-density by \\(1-\alpha)\log 10\\ per
step – and \\\lambda\\ is chosen from it. A Gamma of shape \\0.4\\, for
instance, is sampled through \\X = Y^{1/4}\\.

A density diverging at *both* edges, such as a Beta with both shapes
below one, is beyond any single power: flattening one edge steepens the
other. What does work is a map that behaves like a different power at
each end, \$\$T(u) = \frac{u^{p}}{u^{p} + (1-u)^{q}}, \qquad Y = a +
(b-a)\\T(U),\$\$ monotone from 0 to 1, with a closed-form derivative,
and carrying the exponents \\p\alpha - 1\\ and \\q\beta - 1\\ at the two
ends. Choosing \\p\alpha \> 1\\ and \\q\beta \> 1\\ makes the density of
\\U\\ vanish at both ends, turning a U-shaped density into a
single-peaked one.

Near an edge the variable itself stops being resolvable: at a non-zero
edge the spacing of doubles is absolute, so a slab of \\u\\ collapses
onto the edge and the density evaluated there returns infinity. The
distance from the edge is still exact on the \\u\\ scale, so the power
law is used there instead, calibrated once where the density can be
evaluated. This keeps the mass rather than rejecting it. It cannot place
it any more finely than the arithmetic allows – for a Beta(0.9, 0.1)
some 2.5\\ unit in the last place of 1, and those draws come back equal
to 1 exactly, which is also what
[`rbeta`](https://rdrr.io/r/stats/Beta.html) does.

## Requirements

The kernel must be unimodal and the parameters in `theta` must be
scalars. Densities that diverge at one or at both edges of their support
are handled by the reparameterisations described below.

Heavy tails are not an obstacle: with the default `r = 2` the sampler
handles a Student's t with half a degree of freedom and a Pareto with
infinite mean, and multimodal densities are usually accepted as well.

## See also

[`distrib_rng`](https://statmodels7.github.io/distributions7/reference/distrib_rng.md),
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)

## Examples

``` r
# A distribution defined by its density alone still gets a fast sampler
MyDist <- S7::new_class("MyDist", parent = continuous_distrib)
S7::method(distrib_pdf, MyDist) <- function(distrib, y, theta, log = FALSE) {
  ld <- -log(2 * theta$b) - abs(y - theta$mu) / theta$b
  if (log) ld else exp(ld)
}
d <- MyDist(
  distrib_name = "my laplace", dimension = "univariate", bounds = c(-Inf, Inf),
  params = c("mu", "b"), params_interpretation = c(mu = "location", b = "scale"),
  n_params = 2, params_bounds = list(mu = c(-Inf, Inf), b = c(0, Inf)),
  link_params = list(
    mu = linkfunctions7::identity_link(),
    b = linkfunctions7::log_link()
  )
)
y <- rng_grou(d, 1000, list(mu = 2, b = 1.5))
c(mean = mean(y), var = var(y)) # ~ 2 and ~ 2 * 1.5^2
#>     mean      var 
#> 1.936598 4.797486 
```
