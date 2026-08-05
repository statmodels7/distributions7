# From the Centered Parameters to the Direct Ones

The map \\(\mu, \sigma, \gamma_1) \mapsto (\xi, \omega, \alpha)\\ that
takes the mean, the standard deviation and the skewness to the location,
the scale and the shape.

## Usage

``` r
sn_cp_to_dp(mu, sigma, gamma1, s)
```

## Arguments

- mu, sigma, gamma1:

  The centered parameters, numbers or jets.

- s:

  The sign of \\\gamma_1\\, taken from its plain value.

## Value

A named list with `mu`, `sigma` and `alpha`, the parent's parameters.

## Details

With \\b = \sqrt{2/\pi}\\, \$\$c = \mathrm{sign}(\gamma_1)
\left(\dfrac{2\|\gamma_1\|}{4-\pi}\right)^{1/3}, \qquad \mu_z =
\dfrac{c}{\sqrt{1+c^2}}, \qquad \delta = \dfrac{\mu_z}{b}, \qquad \alpha
= \dfrac{\delta}{\sqrt{1-\delta^2}},\$\$ and then \\\omega =
\sigma/\sqrt{1-\mu_z^2}\\ and \\\xi = \mu - \omega\mu_z\\.

The function is written once and used twice: on plain numbers for the
density, and on jets for the derivatives. The sign is taken from the
plain value before any jet is seeded, which is what a jet cannot do for
itself and what makes this a family rather than a
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md).

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
