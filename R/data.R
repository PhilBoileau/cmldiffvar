#' Toy Population Data
#'
#' A population-level dataset generated according to a simple data-generating
#' process. It provides population-level data on which the functions of this
#' package are tested.
#'
#' @format ## `toy_population_tbl`
#' A [tibble::tibble] with 1,000,000 rows and five columns:
#' \describe{
#' \item{confounder}{Confounding variable of the treatment--outcome
#' relationship} \item{treatment}{A binary treatment assignment indicator}
#' \item{outcome}{A continuous outcome variable}
#' \item{potential_outcome_treament}{A continuous variable representing the
#' potential outcome under treatment assignment}
#' \item{potential_outcome_control}{A continuous variable representing the
#' potential outcome under control assignment} }
#'
"toy_population_tbl"
