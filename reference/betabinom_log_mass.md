# Log-Mass of the Beta-Binomial

The log-mass \\\log\binom{n}{y} + \log B(y+\alpha, n-y+\beta) - \log
B(\alpha, \beta)\\ by whichever of two routes is accurate at the shapes
given.

## Usage

``` r
betabinom_log_mass(y, a, b, n)
```

## Arguments

- y:

  A numeric vector of counts, already known to lie on the support.

- a, b:

  The two shapes.

- n:

  The size.

## Value

A numeric vector of log-probabilities.

## Details

The two beta functions are of magnitude
\\(\alpha+\beta)\log(\alpha+\beta)\\ and their difference is of order
one, so the ordinary route carries an absolute error of \\\varepsilon\\
times that magnitude and is used only while this stays below `1e-8`.
Beyond it the shifts are integers, so each log-gamma difference is an
exact sum of logarithms, \$\$\log\Gamma(\alpha+y) - \log\Gamma(\alpha) =
\sum\_{j=0}^{y-1}\log(\alpha+j),\$\$ and the mass follows from three
such sums without forming any quantity larger than
\\n\log(\alpha+\beta)\\. The sums also give the binomial limit correctly
as the shapes grow at a fixed ratio.
