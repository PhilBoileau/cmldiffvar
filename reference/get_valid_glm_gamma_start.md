# Valid starting value searcher

Automatically look for valid starting values for Gamma GLMs, since Gamma
assumes that \\X\beta \> 0\\.

## Usage

``` r
get_valid_glm_gamma_start(formula, data, max_tries = 100, verbose = FALSE)
```

## Arguments

- formula:

  A `formula` object describing the fit.

- data:

  The data used for the formula.

- max_tries:

  Maximum number of iterations to find valid starting values, default is
  100.

- verbose:

  Flag to determine if the function should print the attempt at which it
  found valid starting values, default is FALSE.

## Value

`start_vals` containing valid starting values.
