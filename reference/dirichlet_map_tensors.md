# Derivative Tensors of a Dirichlet's Shape Vector

\\\alpha = \phi\\\mu(\eta)\\ is bilinear in the concentration and the
mean, so \\\partial^{S}\alpha\\ is \\\phi\\\partial^{S}\mu\\ when \\S\\
names no \\\phi\\, is \\\partial^{S'}\mu\\ when it names one, and is
zero when it names two or more.

## Usage

``` r
dirichlet_map_tensors(s, eta, phi, order)
```

## Arguments

- s:

  A parameters7 `simplex` parameter.

- eta:

  The mean's free vector.

- phi:

  The concentration.

- order:

  The highest order required, 1 to 4.

## Value

A named list of numeric vectors over the coordinates, keyed by sorted
tuple over the composite index set, \\\phi\\ last.
