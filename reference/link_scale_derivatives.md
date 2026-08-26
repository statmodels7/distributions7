# Derivatives on the Link (Real) Scale

The derivative generics of the package
([`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md))
accept a `scale` argument selecting the parametrization the derivatives
are taken with respect to:

- `scale = "parameter"` (default): derivatives with respect to the
  parameters \\\theta\\ on their natural, possibly constrained scale.

- `scale = "link"`: derivatives with respect to the unconstrained linear
  predictors \\\eta_i = g_i(\theta_i)\\, where \\g_i\\ is the link
  stored in `distrib@link_params`. This is the scale on which
  optimization is normally carried out, since \\\eta \in \mathbb{R}^p\\.

## Value

Nothing. This page documents the `scale` argument shared by the
derivative generics named above; the value returned is theirs.

## Details

Write \\h_i = g_i^{-1}\\, so that \\\theta_i = h_i(\eta_i)\\. Because
each parameter carries its own link, the Jacobian
\\\partial\theta/\partial\eta\\ is **diagonal** and the multivariate Faa
di Bruno formula factorizes. For a derivative whose multi-index involves
the distinct parameters \\p_1,\dots,p_r\\ with multiplicities
\\m_1,\dots,m_r\\, \$\$\frac{\partial^k \ell}{\partial \eta\_{p_1}^{m_1}
\cdots \partial \eta\_{p_r}^{m_r}} = \sum\_{j_1=1}^{m_1} \cdots
\sum\_{j_r=1}^{m_r} \ell\_{\\p_1^{j_1} \cdots p_r^{j_r}}
\prod\_{t=1}^{r} B\_{m_t, j_t}\\\left(h\_{p_t}', h\_{p_t}'',
\dots\right)\$\$ where \\\ell\_{\cdot}\\ are the parameter-scale
derivatives and \\B\_{m,j}\\ are the partial (incomplete) Bell
polynomials. The low-order cases are \$\$\frac{\partial \ell}{\partial
\eta_a} = \ell_a h', \qquad \frac{\partial^2 \ell}{\partial \eta_a^2} =
\ell\_{aa} (h')^2 + \ell_a h''\$\$ \$\$\frac{\partial^3 \ell}{\partial
\eta_a^3} = \ell\_{aaa}(h')^3 + 3\ell\_{aa} h' h'' + \ell_a h'''\$\$
\$\$\frac{\partial^4 \ell}{\partial \eta_a^4} = \ell\_{aaaa}(h')^4 +
6\ell\_{aaa}(h')^2 h'' + \ell\_{aa}\left(4h'h''' + 3(h'')^2\right) +
\ell_a h''''\$\$ while mixed derivatives across *different* parameters
simply multiply the corresponding factors (e.g.
\\\partial^2\ell/\partial\eta_a\partial\eta_b = \ell\_{ab} h_a' h_b'\\).

For the **expected** Hessian the first-order term vanishes because the
score has zero expectation, so the information transforms as the
congruence \\\mathrm{diag}(h')\\\mathbb{E}\[H\]\\\mathrm{diag}(h')\\.

The inverse-link derivatives \\h', h'', h''', h''''\\ are obtained from
[`linkfunctions7::linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.html),
so link-scale derivatives are available up to order 4 for every link in
linkfunctions7.

## See also

[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md),
the five generics that take `scale`;
[`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md)
for the polynomials this assembles;
[`linkfunctions7::linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.html),
which supplies the four inverse-link derivatives;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which optimizes on this scale for the reason given above.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 2)
th <- list(mu = 0, sigma = 2)

# The scale carries a log link, so h' = sigma and the first-order case is
# one multiplication.
rbind(parameter = distrib_gradient(d, y, th)$sigma,
      link = distrib_gradient(d, y, th, scale = "link")$sigma)
#>             [,1] [,2] [,3]
#> parameter -0.375 -0.5    0
#> link      -0.750 -1.0    0

# At second order the term in h'' appears, so the link-scale component is
# NOT the parameter-scale one times h' squared.
g <- distrib_gradient(d, y, th)$sigma
h <- distrib_hessian(d, y, th)$sigma_sigma
all.equal(distrib_hessian(d, y, th, scale = "link")$sigma_sigma,
          h * 2^2 + g * 2)
#> [1] TRUE

# For the EXPECTED Hessian the first-order term vanishes, the score having
# mean zero, so the information transforms by congruence.
e <- distrib_expected_hessian(d, y, th)$sigma_sigma
all.equal(distrib_expected_hessian(d, y, th, scale = "link")$sigma_sigma,
          e * 2^2)
#> [1] TRUE

# A parameter on the identity link is unchanged at every order.
all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu,
          distrib_deriv4(d, y, th, scale = "link")$mu_mu_mu_mu)
#> [1] TRUE
```
