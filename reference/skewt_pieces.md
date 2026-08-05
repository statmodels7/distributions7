# The Pieces a Skew t Evaluates From

Assembles the standardized variable and the argument of the tilting
distribution function, together with the six scalar functions every
derivative of the log-density is a combination of.

## Usage

``` r
skewt_pieces(y, mu, sigma, alpha, nu)
```

## Arguments

- y:

  A numeric vector of observations.

- mu, sigma, alpha, nu:

  The parameters.

## Value

A list with `z`, `w`, `c`, `a`, `da`, `e`, `b`, `db`, `q` and `dq`.

## Details

With \\z = (y-\mu)/\sigma\\, \\m = \nu + 1\\ and \\c = \sqrt{m/(\nu +
z^2)}\\, the tilting argument is \\w = \alpha z c\\. The functions
returned are \$\$A = \dfrac{\partial}{\partial z}\log t\_\nu(z) =
-\dfrac{m z}{\nu + z^2}, \qquad A' = -\dfrac{m(\nu - z^2)}{(\nu +
z^2)^2},\$\$ \$\$E = \dfrac{\nu\sqrt{m}}{(\nu + z^2)^{3/2}}, \qquad B =
\dfrac{\partial w}{\partial z} = \alpha E, \qquad B' =
-\dfrac{3\alpha\nu\sqrt{m}\\z}{(\nu + z^2)^{5/2}},\$\$ and \\Q =
t_m(w)/T_m(w)\\ with \\Q' = Q\left\\-\dfrac{(m+1)w}{m + w^2} -
Q\right\\\\, the last from differentiating the quotient.
