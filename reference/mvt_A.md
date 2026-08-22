# The Degrees of Freedom at Large nu, Without the Cancellation

\\A_p(\nu) = \psi(\tfrac{\nu+p}{2}) - \psi(\tfrac{\nu}{2}) -
\tfrac{p}{\nu}\\ and the corresponding first-order combination of the
observed Hessian, each written so that the cancellation between their
terms is performed EXACTLY rather than left to double precision.

## Usage

``` r
mvt_A(nu, p)

mvt_A1(nu)

mvt_T(nu, p)

mvt_T1(nu)

mvt_D(u)
```

## Arguments

- nu:

  The degrees of freedom.

- p:

  The dimension.

## Value

A numeric vector the length of `nu`.

## Details

As \\\nu\\ grows the multivariate t tends to the multivariate gaussian
and every derivative in \\\nu\\ vanishes, so each is a difference of
terms that agree to leading order: \\\psi(\tfrac{\nu+p}{2}) -
\psi(\tfrac{\nu}{2})\\ is \\p/\nu\\ and so is what it is subtracted
from. Measured on the score the family returns, the direct form is wrong
by 4.9e-05 at \\\nu = 10^6\\, by 0.397 at \\10^8\\ and by 838 at
\\10^{10}\\, and it CHANGES SIGN: at \\\nu = 10^9\\ it reads +8.2e-15
where the trend of the values below it gives -2.2e-16.

The repair needs no series, because \\p\\ is an integer dimension and
the shift between the two arguments is therefore a whole number of steps
of the recurrence \\\psi(x+1) = \psi(x) + 1/x\\. For even \\p\\ that
gives sums whose terms all carry ONE SIGN, so nothing cancels:
\$\$A_p(\nu) = -\sum\_{j=0}^{p/2-1} \frac{4j}{\nu(\nu+2j)},\$\$
\$\$\tfrac12\Big\[\psi'(\tfrac{\nu+p}{2}) -
\psi'(\tfrac{\nu}{2})\Big\] + \frac{p}{\nu^2} = \sum\_{j=0}^{p/2-1}
\frac{8j(\nu+j)}{\nu^2(\nu+2j)^2}.\$\$ Both are exactly zero at \\p =
2\\, which is what the two sides of each identity give there and what
the direct forms return as noise at 1e-16.

For odd \\p\\ the shift is a half-integer, so the recurrence carries the
quantity onto the UNIVARIATE \\A_1(\nu)\\, which keeps a series of its
own above a measured crossover – the same expansion distributions7's own
`student_t.cpp` carries, and the one place in this package where that
series exists twice. The two are pinned against each other in the tests.

## See also

[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MvStudentTDistrib.md)
