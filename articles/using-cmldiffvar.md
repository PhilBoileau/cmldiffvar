# Using the cmldiffvar Package

## Background

The average treatment effect (ATE) is a common causal estimand in causal
analyses. Defined as the difference of the mean potential outcomes, it
corresponds to the difference in average outcomes had all study units
been assigned to the treatment groups versus had they all been assigned
to the control group. While the ATE provides an important summary of the
treatment effect, it says nothing of effect heterogeneity. Even when the
ATE indicates a net positive effect of treatment on the average outcome,
some study units may experience better outcomes when assigned to the
control. Decisions informed solely by inference about the ATE might
therefore produce unwelcome results; effect heterogeneity, if it exists,
must be accounted for.

The conditional ATE (CATE), another causal parameter, has gained
popularity as it can uncover and quantify this effect heterogeneity. As
its name suggests, the CATE reports the ATE of study units as a function
of their pre-treatment covariates. However, the CATE can only capture
treatment effect heterogeneity if the considered covariates are
treatment effect modifiers. When the true effect modifiers are unknown
or they are missing from the data due to, for example, resource
constraints or measurement issues, then the CATE will incorrectly
suggest that a treatment effect is homogeneous. Another inferential
strategy not subject to these limitations must supplement inference
about the ATE.

We propose causal inference about the potential outcomes’ variances to
accomplish just that. In particular, we have shown that tests about
contrasts of potential outcomes’ variances, which we refer to
differential variance, can be used to formally test the homogeneous
treatment effect assumption (Boileau et al. In preparation). We briefly
introduce this procedure in the next section.

## Differential Variance Inference

### Complete Data Model

We assume that data are generated according to a parallel group study
design wherein study units are randomly assigned to one of two groups,
possibly based on their pre-treatment covariates. Under a complete data
model — that is, a causal model — we assume that we have access to $n$
independent and identically distributed (i.i.d.) copies of
$X = \left( W,A,Y^{(0)},Y^{(1)} \right)$. Here, $W$ is a random vector
of covariates that includes the treatment–outcome confounders, if any.
$A$ is a random binary variable representing treatment assignment; study
units are said to be assigned to the treatment group when $A = 1$ and to
the control group otherwise. $Y^{(0)}$ and $Y^{(1)}$ are potential
outcomes, continuous random variables representing the outcome under
assignment to the treatment and control groups, respectively.

### Causal Estimands

While we will never have access to copies of $X$, the complete data
model allows us to clearly define the causal estimands of interest. In
particular, we will consider differential variance parameters based on
the potential outcome variances, defined as
$\sigma_{c}^{2}(a) = Var\left\lbrack Y^{(a)} \right\rbrack$. The
potential outcome variance under treatment is given by
$\sigma_{c}^{2}(1)$, that under control is given by $\sigma_{c}^{2}(0)$.
Here, the subscript $c$ indicates that these are parameters of a
complete data model.

#### Absolute Differential Variance

The first differential variance that we consider is the difference of
the potential outcomes’ standard deviations. That is,
$$\psi_{c} = \sqrt{\sigma_{c}^{2}(1)} - \sqrt{\sigma_{c}^{2}(0)}\;.$$
This causal estimand contrasts the variability between groups on the
absolute scale. Positive values indicate that the treatment’s potential
outcomes are more variable than the control’s, negative values indicate
the opposite, and a value of zero indicates that the potential outcomes
have identical variance.

#### Relative Differential Variance

The second differential variance considered here is the ratio of the
potential outcomes’ variances:
$$\lambda_{c} = \frac{\sigma_{c}^{2}(1)}{\sigma_{c}^{2}(0)}\;.$$ This
causal estimand reports the groups variability on the relative scale.
Values larger than one correspond to the potential outcomes under
treatment are more variable, values below zero that the potential
outcomes under control are more variable, and a value of one indicates
that the potential outcomes’ variances are identical.

### Uncovering Heterogeneous Treatment Effects

Now, if the treatment effect were homogeneous, $\psi_{c}$ would
necessarily equal zero and $\lambda_{c}$ would equal 1. This is because,
in a homogeneous treatment effect regime, the potential outcomes under
treatment equal the potential outcomes under control shifted by the ATE.
By the linearity of expectation, the variances of the potential outcomes
are therefore identical.

Conversely, if $\psi_{c} \neq 0$ (or, equivalent, if
$\lambda_{c} \neq 1$), the treatment effect cannot be homogeneous.
Causal inference about $\psi_{c}$ and $\lambda_{c}$ therefore provides a
straightforward test about treatment effect heterogeneity. The formal
theorem and accompanying proof of this result are presented provided by
Boileau et al. (In preparation).

### Observed Data and Identifiability Conditions

The differential variance estimands introduced above are defined in
terms of potential outcomes which are never available in practice.
Luckily, we can rely on standard — but unverifiable — causal
identifiability assumptions to perform inference about these causal
estimands from observable data.

Under an observed data model, we have access to $n$ i.i.d. copies of
random observations $O = (W,A,Y)$, where $W$ and $A$ are as defined in
$X$. Here, $Y$ corresponds to a random variable of the observation’s
continuous outcome. Then, under the assumption of full exchangeability,
positivity, and consistency, we obtain for $a = 0$ or $a = 1$ the
following results (Boileau et al. In preparation): $$\begin{aligned}
{\sigma_{c}^{2}(a)} & {= {\mathbb{E}}\left\lbrack {\mathbb{E}}\left\lbrack Y^{2}|W,A = a \right\rbrack \right\rbrack - {\mathbb{E}}\left\lbrack {\mathbb{E}}\left\lbrack Y|W,A = a \right\rbrack \right\rbrack^{2}} \\
 & {= \sigma^{2}(a)\;.}
\end{aligned}$$ The assumption of full exchangeability assumes that $W$
contains all confounders of the treatment–outcome relationship. The
positive assumption requires that study units have a non-zero
probability of being assigned to either group. Finally, the consistency
assumptions states that the observed outcomes correspond to the
potential outcome associated with the study units’ group assignments. Of
note, these identifiability conditions are satisfied by construction in
randomized controlled trials. For an in-depth discussion of these
assumptions, see Hernan and Robins (2025).

When these conditions are satisfied, we can in turn define parameters in
the observed data model that are equivalent to the differential variance
estimands previously defined in the complete data model. That is.
$$\begin{aligned}
\psi_{c} & {= \sqrt{\sigma^{2}(1)} - \sqrt{\sigma^{2}(0)}} \\
 & {= \psi}
\end{aligned}$$ and $$\begin{aligned}
\lambda_{c} & {= \frac{\sigma^{2}(1)}{\sigma^{2}(0)}} \\
 & {= \lambda\;.}
\end{aligned}$$

We can therefore causal perform inference about these differential
variance parameters without requiring access to the potential outcomes,
assuming the causal identifiability conditions are satisfied. The
`cmldiffvar` package implements such procedures, which we briefly
describe next.

## `cmldiffvar` in Action

The `cmldiffvar` package implements two types of estimators for the
absolute and relative differential variance: unadjusted and adjusted
estimators. The former is appropriate for randomized experiments in
which there are no confounders of the treatment–outcome relationship.
The latter should be used in analyses of observational study data to
account for confounding variables. Adjusted estimators can also be used
in randomized experiments to adjust for pre-treatment covariates that
explain variability in the outcome variable, potentially improving upon
the precision of the unadjusted estimator.

### Unadjusted Estimators

The unadjusted estimators of the absolute and relative differential
variance are defined as the difference of the groups’ outcomes’ standard
deviations and the ratio of the groups’ outcomes’ variances,
respectively. That is,
$$\widehat{\psi} = \sqrt{{\widehat{\sigma}}^{2}(1)} - \sqrt{{\widehat{\sigma}}^{2}(0)}$$
and
$$\widehat{\lambda} = \frac{{\widehat{\sigma}}^{2}(1)}{{\widehat{\sigma}}^{2}(0)}\;,$$
where
$${\widehat{\sigma}}^{2}(a) = \widehat{\text{Var}}\left\lbrack Y|A = a \right\rbrack\;,$$
the group-specific sample variance.

The unadjusted absolute and relative differential variance estimators
are implemented in the
[`unadjdiffvar()`](https://pboileau.ca/cmldiffvar/reference/unadjdiffvar.md)
function. The absolute differential variance estimator is called by
default; the relative differential variance estimator is called by
setting the `estimand_type` function to `"relative"`.

These estimators’ asymptotic sampling distributions are normal with
means equal to their respective target estimand and with variance given
by the variances of their respective influence function. Their influence
functions can therefore be used to construct Wald-type 95% confidence
intervals, which in turn permit hypothesis testing. In particular, these
estimators can be used to perform hypothesis tests about the homogeneous
treatment effect assumption. Under the null hypothesis, $\psi = 0$ and
$\lambda = 1$. Under the alternative, $\psi \neq 0$ and
$\lambda \neq 1$.

#### Example

We estimate the relative differential variance from a simulated
randomized controlled trial using the
[`unadjdiffvar()`](https://pboileau.ca/cmldiffvar/reference/unadjdiffvar.md)
function. We must specify the estimand using the `estimand_type`
argument, and provide the variable names of the treatment, outcome, and
propensity score variable names in the sample data set. If the
propensity score variable name is not provided to
[`unadjdiffvar()`](https://pboileau.ca/cmldiffvar/reference/unadjdiffvar.md),
they will be estimated from the provided data set as the proportion of
study units assigned to the treatment group.

The variances of the potential outcomes under the treatment and control
group assignments are $1.44$ and $1$, respectively. The relative
differential variance parameter therefore has a value of
$\lambda_{c} = 1.44$, indicating that the treatment effect is
heterogeneous.

``` r
# load the required libraries
library(cmldiffvar)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# set the seed for reproducibility
set.seed(3976468)

# simulate randomized controlled trial population
n_pop <- 100000
propensity_score <- 0.5
propensity_score_vec <- rep(propensity_score, n_pop)
treatment_vec <- rbinom(n_pop, 1, propensity_score)
outcome_treatment_vec <- rnorm(n = n_pop, mean = 3, sd = 1.2)
outcome_control_vec <- rnorm(n = n_pop, mean = 1, sd = 1)
outcome_vec <- treatment_vec * outcome_treatment_vec +
  (1 - treatment_vec) * outcome_control_vec
population_tbl <- tibble(
  propensity_score = propensity_score_vec,
  treatment = treatment_vec,
  outcome = outcome_vec
)

# randomly sample from the population
sample_tbl <- slice_sample(population_tbl, n = 500)

# estimate the relative differential variance and test the homogeneous treatment
# assumption
rel_diff_var_est <- sample_tbl |>
  unadjdiffvar(
    estimand_type = "relative",
    treatment_var_name = "treatment",
    propensity_score_var_name = "propensity_score",
    outcome_var_name = "outcome"
  )

rel_diff_var_est
#> # A tibble: 1 × 6
#>   estimand                       estimate    se ci_low ci_high p_value
#>   <chr>                             <dbl> <dbl>  <dbl>   <dbl>   <dbl>
#> 1 relative differential variance     1.51 0.187   1.14    1.87 0.00676
```

The unadjusted estimator produces an estimate of 1.507, suggesting that
the potential outcomes under the treatment condition are more variable
than the potential outcomes under control. The 95% confidence interval
does not cover $1$, providing evidence of a statistically significant
heterogeneous treatment effect. Indeed, the p-value of the associated
test is 0.007 such that the homogeneous treatment effect hypothesis is
rejected.

### Adjusted Estimators

The `cmldiffvar` package implements two adjusted estimators for each
differential variance parameter: a one-step estimator (Bickel et al.
1993; Pfanzagl and Wefelmeyer 1985) and a targeted maximum likelihood
estimator (TMLE) (van der Laan and Rubin 2006; van der Laan and Rose
2011, 2018). Note that the one-step estimator is equivalent to an
estimating equation (alternatively known as a double/debiased machine
learning estimator) (van der Laan and Robins 2003; Chernozhukov et al.
2017). These estimators rely on the estimation of three nuisance
parameters:

1.  the propensity score: $g(W) = P\left( A = 1|W \right)$,
2.  the outcome regression:
    $\bar{Q}(W,A) = E\left\lbrack Y|W,A \right\rbrack$, and
3.  the squared outcome regression:
    ${\bar{Q}}^{2}(W,A) = E\left\lbrack Y^{2}|W,A \right\rbrack$ .

These adjusted estimators are doubly robust. That is, if $g(W)$ or both
$\bar{Q}(W,A)$ and ${\bar{Q}}^{2}(W,A)$ are consistently estimated, then
the adjusted estimators are consistent as well. Further, if the nuisance
estimators of $g(W)$, $\bar{Q}(W,A)$, and ${\bar{Q}}^{2}(W,A)$ are
consistent and converge at a rate of $n^{1/4}$ under the $L_{2}$ norm,
then these estimators are regular and asymptotically linearly. This
implies that these differential variance parameter estimators’ sampling
distributions are asymptotically normal about the true differential
variance parameters with variance given by the parameters’ efficient
influence functions. As with the unadjusted estimators, 95% Wald-type
confidence intervals can in turn be constructed, and hypotheses about
treatment effect heterogeneity can be performed. Additional details are
provided in Boileau et al. (In preparation).

Now, this package relies on Super Learners to flexibly estimate $g(W)$,
$\bar{Q}(W,A)$, and ${\bar{Q}}^{2}(W,A)$(van der Laan, Polley, and
Hubbard 2007). In particular, the `cmldiffvar` package relies on the
[`SuperLearner`](https://cran.r-project.org/web/packages/SuperLearner/index.html)
R package implementation of this ensemble learning framework. By
default, a generalized linear model, multivariate adaptive regression
splines (Friedman 1991), and a Random Forest (Breiman 2001) are used as
base learners in the Super Learner estimator of the propensity score and
the outcome regression. A generalized quadratic linear model and a
Random Forest make up the base learners of the Super Learner estimator
for the squared outcome regression. Of course, an alternative set of
base learners used for each estimator can be specified in
[`cmldiffvar()`](https://pboileau.ca/cmldiffvar/reference/cmldiffvar.md)
function.

Finally, we note that these adjusted estimators are guaranteed to be
regular and asymptotically linear when applied to data generated in a
randomized controlled trial. The known propensity scores, corresponding
to study units arm randomization probabilities, can be provided to the
[`cmldiffvar()`](https://pboileau.ca/cmldiffvar/reference/cmldiffvar.md)
function using its `propensity_score_var_name` argument.

#### Example

We estimate the absolute differential variance from a random of a
simulated observational study population provided with the `cmldiffvar`
package using the TMLE implemented in the
[`cmldiffvar()`](https://pboileau.ca/cmldiffvar/reference/cmldiffvar.md)
function. Similar to the unadjusted estimator, we must provide the names
of the treatment and outcome variables in the sample data set. The
treatment–outcome relationship is confounded by a single variable called
`"confounder"`, which we adjust for using the
`propensity_score_adj_var_names` and `cond_exp_outcome_adj_var_names`
arguments. Note that the adjustment sets for the propensity score and
conditional (squared) outcome regressions are specified separately to
allow for increased flexibility when estimating nuisance parameters.
Finally, we rely on the default Super Learners to estimate the nuisance
parameters.

The variance of the potential outcomes under treatment and control are
$10$ and $2$, respectively. The absolute differential variance estimand
therefore has a value of $1.75$.

``` r
# load the required libraries
library(cmldiffvar)
library(dplyr)
library(SuperLearner)
#> Loading required package: nnls
#> Loading required package: gam
#> Loading required package: splines
#> Loading required package: foreach
#> Loaded gam 1.22-6
#> Super Learner
#> Version: 2.0-29
#> Package created on 2024-02-06

# set the seed for reproducibility
set.seed(86523)

# randomly sample from the provided observational study population
sample_tbl <- slice_sample(toy_population_tbl, n = 500)

# estimate the relative differential variance and test the homogeneous treatment
# assumption
abs_diff_var_est <- sample_tbl |>
  cmldiffvar(
    estimator_type = "tmle",
    treatment_var_name = "treatment",
    outcome_var_name = "outcome",
    propensity_score_adj_var_names = "confounder",
    cond_exp_outcome_adj_var_names = "confounder" 
  )
#> Loading required namespace: earth
#> Loading required namespace: ranger

abs_diff_var_est
#> # A tibble: 1 × 6
#>   estimand                       estimate    se ci_low ci_high  p_value
#>   <chr>                             <dbl> <dbl>  <dbl>   <dbl>    <dbl>
#> 1 absolute differential variance     1.62 0.152   1.32    1.92 1.28e-26
```

The TMLE produces an estimate of 1.622, indicating that the standard
deviation of the potential outcomes under treatment is greater than that
of the potential outcomes under control. Further, the 95% confidence
interval and the p-value of 0 suggest that we reject the null hypothesis
of a homogeneous treatment effect. Near identical results are obtained
when using a one-step estimator instead, which can be specified by
setting the `estimator_type` argument to `"one-step"`.

## Session Information

    #> R version 4.5.2 (2025-10-31)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.3 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] splines   stats     graphics  grDevices utils     datasets  methods  
    #> [8] base     
    #> 
    #> other attached packages:
    #> [1] SuperLearner_2.0-29 gam_1.22-6          foreach_1.5.2      
    #> [4] nnls_1.6            dplyr_1.1.4         cmldiffvar_0.0.1   
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] Matrix_1.7-4      jsonlite_2.0.0    ranger_0.17.0     compiler_4.5.2   
    #>  [5] Rcpp_1.1.0        plotrix_3.8-4     tidyselect_1.2.1  jquerylib_0.1.4  
    #>  [9] systemfonts_1.3.1 textshaping_1.0.4 yaml_2.3.10       fastmap_1.2.0    
    #> [13] lattice_0.22-7    R6_2.6.1          generics_0.1.4    Formula_1.2-5    
    #> [17] knitr_1.50        iterators_1.0.14  backports_1.5.0   checkmate_2.3.3  
    #> [21] tibble_3.3.0      desc_1.4.3        bslib_0.9.0       pillar_1.11.1    
    #> [25] rlang_1.1.6       utf8_1.2.6        cachem_1.1.0      xfun_0.54        
    #> [29] fs_1.6.6          sass_0.4.10       earth_5.3.4       cli_3.6.5        
    #> [33] pkgdown_2.2.0     withr_3.0.2       magrittr_2.0.4    grid_4.5.2       
    #> [37] digest_0.6.37     plotmo_3.6.4      lifecycle_1.0.4   vctrs_0.6.5      
    #> [41] evaluate_1.0.5    glue_1.8.0        codetools_0.2-20  ragg_1.5.0       
    #> [45] rmarkdown_2.30    tools_4.5.2       pkgconfig_2.0.3   htmltools_0.5.8.1

## References

Bickel, Peter J., Chris A. J. Klaassen, Ya’acov Ritov, and Jon A.
Wellner. 1993. *Efficient and Adaptive Estimation for Semiparametric
Models*. Springer New York.

Boileau, Philippe A, Hani Zaki, Gabriele Lileikyte, Niklas Nielsen,
Patrick R Lawler, and Mireille E Schnitzer. In preparation.
“Assumption-Lean Differential Variance Inference for Heterogeneous
Treatment Effect Detection.”

Breiman, Leo. 2001. “Random Forests.” *Machine Learning* 45 (1): 5–32.
<https://doi.org/10.1023/A:1010933404324>.

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo,
Christian Hansen, and Whitney Newey. 2017. “Double/Debiased/Neyman
Machine Learning of Treatment Effects.” *American Economic Review* 107
(5): 261–65. <https://doi.org/10.1257/aer.p20171038>.

Friedman, Jerome H. 1991. “Multivariate Adaptive Regression Splines.”
*The Annals of Statistics* 19 (1): 1–67.
<https://doi.org/10.1214/aos/1176347963>.

Hernan, Miguel A., and James M. Robins. 2025. *Causal Inference: What
If*. CRC Press.

Pfanzagl, Johann, and Wolfgang Wefelmeyer. 1985. “Contributions to a
General Asymptotic Statistical Theory.” *Statistics & Risk Modeling* 3
(3-4): 379–88. <https://doi.org/10.1524/strm.1985.3.34.379>.

van der Laan, Mark J., Eric C. Polley, and Alan E. Hubbard. 2007. “Super
Learner.” *Statistical Applications in Genetics and Molecular Biology* 6
(1). <https://doi.org/10.2202/1544-6115.1309>.

van der Laan, Mark J., and James M. Robins. 2003. *Unified Methods for
Censored Longitudinal Data and Causality*. Springer Series in
Statistics. New York, NY: Springer.
<https://doi.org/10.1007/978-0-387-21700-0>.

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
