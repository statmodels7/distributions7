# Assemble a Multivariate Student t's Higher Derivatives

Splits a component into its mean-and-matrix part and its \\\nu\\ part.
The log-determinant and the \\\nu\\ constant each contribute to one kind
of component only; the remaining term is a univariate chain rule in
\\q\\ whose outer derivatives, \\(-1)^{m-1}(m-1)!/(\nu+q)^{m}\\, are
then differentiated in \\\nu\\ by the Leibniz rule against the prefactor
\\(\nu+p)/2\\, which is linear.

## Usage

``` r
mvt_higher(distrib, y, theta, order)
```

## Arguments

- distrib:

  A `MvStudentTDistrib` object.

- y:

  An \\n \times p\\ matrix of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 to 4.

## Value

A named list of numeric vectors, one per component.
