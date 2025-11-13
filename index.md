# R/`cmldiffvar`

> Causal Machine Learning Methods for Differential Variance Inference

**Authors:** [Philippe Boileau](https://pboileau.ca), [Hani
Zaki](https://ca.linkedin.com/in/hani-zaki), [Mireille
Schnizter](https://www.mireilleschnitzer.com)

------------------------------------------------------------------------

## What’s `cmldiffvar`?

`cmldiffvar` implements causal machine learning methods for differential
variance inference. These methods rely on semiparametric efficiency
theory and flexible machine learning methods — namely, Super Learner
ensembles — to avoid the need for convenience assumptions about
data-generating processes (van der Laan and Rose 2011; van der Laan,
Polley, and Hubbard 2007). Hypothesis tests about differential variance
can uncover heterogeneous treatment effects, even when the effect
modifiers are excluded from the data. Details on the methodology are
provided in Boileau et al. (In preparation).

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
difference of the potential outcomes’ standard deviations, on a random
sample of the `toy_population_tbl` data included with the `cmldiffvar`
package. This dataset represents an observational study in which the
treatment variable is binary, the outcome is continuous, and a single
confounder was measured. The true absolute differential variance in this
population is $2$. Because the absolute differential variance is
non-zero, the treatment effect is heterogeneous.

We use a targeted maximum likelihood estimator (van der Laan and Rubin
2006; van der Laan and Rose 2011, 2018), the
[`cmldiffvar()`](https://pboileau.ca/cmldiffvar/reference/cmldiffvar.md)
function’s default estimator, to infer the differential variance of this
population. The function outputs a point estimate and a $95\%$
confidence interval by default. A p-value corresponding to a test of
whether the differential variance is significantly different from zero
is also provided. This is equivalent to testing whether the treatment
effect is homogeneous.

``` r
# load the required packages
library(cmldiffvar)
library(dplyr)
library(SuperLearner)

# set the seed for reproducibility
set.seed(510)

# random sample from population data
sample_tbl <- slice_sample(toy_population_tbl, n = 250)

# estimate absolute differential variance
dif_var_result_tbl <- sample_tbl |>
  cmldiffvar(
    propensity_score_adj_var_names = "confounder",
    cond_exp_outcome_adj_var_names = "confounder",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome"
  )
```

| estimand                       | estimate |   se | ci_low | ci_high | p_value |
|:-------------------------------|---------:|-----:|-------:|--------:|--------:|
| absolute differential variance |     2.11 | 0.29 |   1.53 |    2.69 |    0.00 |

The absolute differential variance point estimate is near the ground
truth. Additionally, the test correctly rejects the null hypothesis of a
homogeneous treatment effect at the $5\%$ significance level.

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

``` R
@unpublished{boileau2025,
 author = {Philippe A Boileau and Hani Zaki and Gabriele Lileikyte and Niklas
           Nielsen and Patrick R Lawler and Mireille E Schnitzer},
 title = {Assumption-Lean Differential Variance Inference for Heterogeneous
          Treatment Effect Detection},
 year = {In preparation}
}
```

------------------------------------------------------------------------

## Licence

© 2025 [Philippe Boileau](https://pboileau.ca)

The contents of this repository are distributed under the MIT license.
See file
[`LICENSE.md`](https://github.com/PhilBoileau/cmldiffvar/blob/main/LICENSE.md)
for details.

------------------------------------------------------------------------

## References

Boileau, Philippe A, Hani Zaki, Gabriele Lileikyte, Niklas Nielsen,
Patrick R Lawler, and Mireille E Schnitzer. In preparation.
“Assumption-Lean Differential Variance Inference for Heterogeneous
Treatment Effect Detection.”

van der Laan, Mark J., Eric C. Polley, and Alan E. Hubbard. 2007. “Super
Learner.” *Statistical Applications in Genetics and Molecular Biology* 6
(1). <https://doi.org/10.2202/1544-6115.1309>.

van der Laan, Mark J., and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Series in
Statistics. New York, NY: Springer.
<https://doi.org/10.1007/978-1-4419-9782-1>.

———. 2018. *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Series in Statistics. Cham:
Springer International Publishing.
<https://doi.org/10.1007/978-3-319-65304-4>.

van der Laan, Mark J., and Daniel Rubin. 2006. “Targeted Maximum
Likelihood Learning.” *The International Journal of Biostatistics* 2
(1). <https://doi.org/10.2202/1557-4679.1043>.
