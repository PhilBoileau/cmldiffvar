#' Toy Population Data
#'
#' A population-level dataset generated according to a simple data-generating
#' process. It provides population-level data on which the functions of this
#' package are tested.
#'
#' @format ## `toy_population_tbl` A [tibble] with 100,000 rows and five columns:
#' \describe{
#' \item{confounder}{Confounding variable of the treatment--outcome
#' relationship}
#' \item{treatment}{A binary treatment assignment indicator}
#' \item{propensity_score}{A numeric variable corresponding to the propensity
#' score}
#' \item{outcome}{A continuous outcome variable}
#' \item{potential_outcome_treatment}{A continuous variable representing the
#' potential outcome under treatment assignment}
#' \item{potential_outcome_control}{A continuous variable representing the
#' potential outcome under control assignment} }
#'
"toy_population_tbl"

#' Toy Population Data with multiple confounders
#'
#' A population-level dataset generated according to a simple data-generating
#' process which includes three confounder variables.
#'
#' @format ## `multiple_confounders_toy_population_tbl`
#' A [tibble] with 100,000 rows and seven columns:
#' \describe{
#' \item{confounder_1}{First confounding variable of the treatment--outcome
#' relationship}
#' \item{confounder_2}{Second confounding variable of the treatment--outcome
#' relationship}
#' \item{confounder_3}{Third confounding variable of the treatment--outcome
#' relationship}
#' \item{treatment}{A binary treatment assignment indicator}
#' \item{outcome}{A continuous outcome variable}
#' \item{potential_outcome_treatment}{A continuous variable representing the
#' potential outcome under treatment assignment}
#' \item{potential_outcome_control}{A continuous variable representing the
#' potential outcome under control assignment} }
#'
"multiple_confounders_toy_population_tbl"

