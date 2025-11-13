# Toy Population Data

A population-level dataset generated according to a simple
data-generating process. It provides population-level data on which the
functions of this package are tested.

## Usage

``` r
toy_population_tbl
```

## Format

`toy_population_tbl` A
[tibble](https://tibble.tidyverse.org/reference/tibble.html) with
100,000 rows and five columns:

- `confounder`: Confounding variable of the treatment–outcome
  relationship

- `treatment`: A binary treatment assignment indicator

- `propensity_score`: A numeric variable corresponding to the propensity
  score `outcome`: A continuous outcome variable

- `potential_outcome_treatment`: A continuous variable representing the
  potential outcome under treatment assignment

- `potential_outcome_control`: A continuous variable representing the
  potential outcome under control assignment
