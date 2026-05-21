#' Melt Wide Time-to-Event Data
#'
#' @description `melt_tte_data_fun()` turns a wide-format
#' [`tibble`][tibble::tibble] into a long-format [`tibble`][tibble::tibble].
#'
#'
#' @param wide_data_tbl A wide format [`tibble`][tibble::tibble].
#' @param baseline_var_names A `character` vector of column names corresponding
#'   to baseline covariates.
#' @param treatment_var_name A `character` corresponding to the exposure
#'   variable.
#' @param outcome_var_name A `character` corresponding to the outcome variable.
#' @param censoring_var_name A `character` indicating the right censoring
#'   indicator variable.
#' @param time_cutoff A `numeric` representing the time point at which to
#'   evaluate the time-to-event parameter.
#'
#' @returns A longitudinal time-to-event [`tibble`][tibble::tibble]. That is,
#'   observations are repeated at all time points between study entry and the
#'   `time_cutoff` argument.
#'
#' @keywords internal
#'
melt_tte_data_fun <- function(
  wide_data_tbl,
  baseline_var_names,
  treatment_var_name,
  outcome_var_name,
  censoring_var_name,
  time_cutoff
) {

  # remove unnecessary variables
  wide_data_tbl <- wide_data_tbl |>
    dplyr::select(
      dplyr::all_of(
        c(baseline_var_names, treatment_var_name, outcome_var_name,
          censoring_var_name)
      )
    )

  # get the unique times in the dataset up to the cutoff time
  unique_times <- wide_data_tbl |>
    dplyr::filter(!!rlang::sym(outcome_var_name) <= time_cutoff) |>
    dplyr::pull(outcome_var_name) |>
    unique() |>
    sort()

  # lengthen each observation
  obs_ls <- lapply(
    seq_len(nrow(wide_data_tbl)),
    function(obs_idx) {

      # extract the observation and repeat it for each time point
      long_obs_tbl <- wide_data_tbl |>
        dplyr::slice(rep(obs_idx, times = length(unique_times)))

      # label each entry by the time, create the longitudinal indicators, add an
      # id, and indicate which times have an event associated with them for use
      # in nuisance estimation
      long_obs_tbl <- long_obs_tbl |>
        dplyr::mutate(
          cmldiffvar_long_time = unique_times,
          R = dplyr::if_else(cmldiffvar_long_time == .data[[outcome_var_name]] &
                               .data[[censoring_var_name]] == 1, 1, 0),
          L = dplyr::if_else(cmldiffvar_long_time == .data[[outcome_var_name]] &
                               .data[[censoring_var_name]] == 0, 1, 0),
          I = dplyr::if_else(
            cmldiffvar_long_time <= .data[[outcome_var_name]], 1, 0),
          J = dplyr::if_else(I == 1 & L == 0, 1, 0),
          cmldiffvar_id = obs_idx
        )

      # retain only the necessary variables
      long_obs_tbl |>
        dplyr::select(
          dplyr::all_of(
            c("cmldiffvar_id", "cmldiffvar_long_time", baseline_var_names,
              treatment_var_name, "R", "L", "I", "J")
          )
        )
    }
  )

  # collapse the observation list into a tibble
  dplyr::bind_rows(obs_ls)

}
