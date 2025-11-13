# Toy Population Data with multiple confounders

A population-level dataset generated according to a simple
data-generating process which includes three confounder variables.

## Usage

``` r
multiple_confounders_toy_population_tbl
```

## Format

`multiple_confounders_toy_population_tbl` A
[tibble](https://tibble.tidyverse.org/reference/tibble.html) with
100,000 rows and seven columns:

- `confounder_1`: First confounding variable of the treatment–outcome
  relationship

- `confounder_2`: Second confounding variable of the treatment–outcome
  relationship

- `confounder_3`: Third confounding variable of the treatment–outcome
  relationship

- `treatment`: A binary treatment assignment indicator

- `outcome`: A continuous outcome variable potential_outcome_treatment:
  A continuous variable representing the potential outcome under
  treatment assignment

- `potential_outcome_control`: A continuous variable representing the
  potential outcome under control assignment
