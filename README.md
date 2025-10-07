
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`cmldiffvar`

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/license/mit/)
[![R-CMD-check](https://github.com/PhilBoileau/cmldiffvar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PhilBoileau/cmldiffvar/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/PhilBoileau/cmldiffvar/graph/badge.svg)](https://app.codecov.io/gh/PhilBoileau/cmldiffvar)
<!-- badges: end -->

> Causal Machine Learning Methods for Differential Variance Inference

**Authors:** [Philippe Boileau](https://pboileau.ca), [Hani
Zaki](link-to-profile), [Mireille
Schnizter](https://www.mireilleschnitzer.com)

------------------------------------------------------------------------

## What’s `cmldiffvar`?

`cmldiffvar` implements causal machine learning methods for differential
variance inference. These methods rely on semiparametric efficiency
theory and flexible machine learning methods — namely, Super Learner
ensembles — to avoid the need for convenience assumptions about
data-generating processes (van der Laan and Rose 2011; van der Laan,
Polley, and Hubbard 2007).

------------------------------------------------------------------------

## Installation

The *development version* of the package may be installed from GitHub
using [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("PhilBoileau/cmldiffvar")
```

------------------------------------------------------------------------

## Example

We estimate the absolute differential variance, defined as the
difference of the potential outcome variances, on a random sample of the
`toy_population_tbl` data included with the `cmldiffvar` package. We use
a targeted maximum likelihood estimator that allows adjustment for
confounding variables. The true absolute differential variance in this
population is $2$, indicating that there is treatment effect
heterogeneity.

``` r
# load the required packages
library(cmldiffvar)
library(dplyr)
library(SuperLearner)
library(knitr)

# set the seed for reproducibility
set.seed(510)

# random sample from population data
sample_tbl <- slice_sample(toy_population_tbl, n = 500)

# estimate absolute differential variance
dif_var_result_tbl <- sample_tbl |>
  cmldiffvar(
    propensity_score_adj_var_names = "confounder",
    cond_exp_outcome_adj_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome"
  )

# print formatted table
dif_var_result_tbl |> kable(digits = 2, format.args = list(nsmall = 2))
```

| estimand                       | estimate |   se | ci_low | ci_high | p_value |
|:-------------------------------|---------:|-----:|-------:|--------:|--------:|
| absolute differential variance |     2.09 | 0.19 |   1.72 |    2.47 |    0.00 |

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/PhilBoileau/cmldiffvar/issues).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/PhilBoileau/cmldiffvar/blob/main/.github/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

Please cite the following paper when using the `cmldiffvar` R software
package.

    @article{boileau-cmldiffvar,
      author = {Philippe boileau and Hani Zaki and Mireille Schnizter},
      journal = {arXiv preprint},
      title = {Causal Machine Learning Methods for Differential Varance Inference},
      url = {NA},
      year = {2025}
    }

------------------------------------------------------------------------

## Licence

© 2025 [Philippe Boileau](https://pboileau.ca)

The contents of this repository are distributed under the MIT license.
See file
[`LICENSE.md`](https://github.com/PhilBoileau/cmldiffvar/blob/main/LICENSE.md)
for details.

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-laanSuperLearner2007" class="csl-entry">

van der Laan, Mark J., Eric C. Polley, and Alan E. Hubbard. 2007. “Super
Learner.” *Statistical Applications in Genetics and Molecular Biology* 6
(1). <https://doi.org/10.2202/1544-6115.1309>.

</div>

<div id="ref-vanderlaanTargetedLearningCausal2011" class="csl-entry">

van der Laan, Mark J., and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Series in
Statistics. New York, NY: Springer.
<https://doi.org/10.1007/978-1-4419-9782-1>.

</div>

</div>
