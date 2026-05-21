library(dplyr)

# function for randomly generating a basic TTE outcome dataset for testing
generate_test_data <- function(
    n_obs = 200,
    null_marginal = FALSE
) {

  # baseline covariate
  w <- rbinom(n = n_obs, size = 1, prob = 0.5)

  # exposure
  prop_score <- plogis(-0.2 + 0.7 * w)
  a <- rbinom(n_obs, 1, prop_score)

  # define hazard functions
  max_time <- 100
  cond_cens_hazard <- function(time, exposure, w) {
    if (exposure == 1) {
      0.005
    } else {
      0.008
    }
  }
  if (null_marginal) {
    cond_surv_hazard <- function(time, exposure, w) {
      (time < max_time) / (75 + exp(2 * w)) + (time == max_time)
    }
  } else {
    cond_surv_hazard <- function(time, exposure, w_1, w_2, w_3) {
      (time < max_time) / (75 + exp(2 * exposure + 2 * w)) + (time == max_time)
    }
  }

  # generate the failure events for t = 1 to max_time
  failure_time_sim <- function(exposure) {
    sapply(
      seq_len(n_obs),
      function(obs) {
        failure_time <- NA
        for (t in 1:max_time) {
          prob <- cond_surv_hazard(t, exposure, w[obs])
          status <- rbinom(1, 1, prob)
          if (status == 1) {
            failure_time <- t
            break
          }
        }
        return(failure_time)
      }
    )
  }
  failure_time_1 <- failure_time_sim(1)
  failure_time_0 <- failure_time_sim(0)


  # generate the censoring events for t = 1 to max_time
  censor_time_sim <- function(exposure) {
    sapply(
      seq_len(n_obs),
      function(obs) {
        censor_time <- NA
        for (t in 1:max_time) {
          prob <- cond_cens_hazard(t, exposure, w[obs])
          status <- rbinom(1, 1, prob)
          if (status == 1) {
            censor_time <- t
            break
          }
        }
        if (is.na(censor_time)) censor_time <- max_time + 1
        return(censor_time)
      }
    )
  }
  censor_time_1 <- censor_time_sim(1)
  censor_time_0 <- censor_time_sim(0)

  # compile the failure and censoring times
  failure_time <- sapply(
    seq_len(n_obs),
    function(obs) {
      if (a[obs] == 1) failure_time_1[obs] else failure_time_0[obs]
    }
  )
  censor_time <- sapply(
    seq_len(n_obs),
    function(obs) {
      if (a[obs] == 1) censor_time_1[obs] else censor_time_0[obs]
    }
  )

  # assess the observed time-to-event and censoring indicator
  time <- sapply(
    seq_len(n_obs),
    function(obs) {
      if (censor_time[obs] < failure_time[obs]) {
        censor_time[obs]
      } else {
        failure_time[obs]
      }
    }
  )
  censoring <- sapply(
    seq_len(n_obs),
    function(obs) if (time[obs] == censor_time[obs]) 1 else 0
  )

  # assemble the data.table
  wide_tbl <- dplyr::tibble(
    w = w,
    prop_score = prop_score,
    a = a,
    time = time,
    censoring = censoring,
    potential_time_0 = failure_time_0,
    potential_time_1 = failure_time_1
  )

  return(wide_tbl)
}

# approximate relevant parameter values at time_cutoff = 50
# set.seed(510)
# pop_sim_tbl <- generate_test_data(n_obs = 100000, null_marginal = TRUE)
# time_cutoff <- 50
# pop_params_tbl <- pop_sim_tbl |>
#   mutate(
#     trunc_potential_time_0 = if_else(
#       potential_time_0 > time_cutoff, time_cutoff, potential_time_0
#     ),
#     trunc_potential_time_1 = if_else(
#       potential_time_1 > time_cutoff, time_cutoff, potential_time_1
#     )
#   ) |>
#   summarise(
#     rmst_0 = mean(trunc_potential_time_0),
#     rmst_1 = mean(trunc_potential_time_1),
#     var_0 = var(trunc_potential_time_0),
#     var_1 = var(trunc_potential_time_1),
#     rmst_diff = rmst_1 - rmst_0,
#     diff_var_ratio = var_1 / var_0
#   )

# recorded values
# $ rmst_0         <dbl> 37.17233
# $ rmst_1         <dbl> 37.23047
# $ var_0          <dbl> 279.903
# $ var_1          <dbl> 278.9438
# $ rmst_diff      <dbl> 0.05814
# $ diff_var_ratio <dbl> 0.9965733

# compute the ratios of truncated survival times
# diff_var_ratio_vec <- sapply(
#   seq(from = 0, to = 10),
#   function(trunc_value) {
#
#     # compute the additive diff var ratio under the homogeneous additive
#     # treatment effect assumption
#     add_diff_var_ratio <- pop_sim_tbl |>
#       mutate(
#         trunc_trunc_potential_time_0 = if_else(
#           potential_time_0 > (time_cutoff - trunc_value),
#           (time_cutoff - trunc_value), potential_time_0
#         ),
#         trunc_potential_time_1 = if_else(
#           potential_time_1 > time_cutoff, time_cutoff, potential_time_1
#         )
#       ) |>
#       summarise(
#         var_0 = var(trunc_trunc_potential_time_0),
#         var_1 = var(trunc_potential_time_1),
#         add_diff_var_ratio = var_1 / var_0
#       ) |>
#       pull(add_diff_var_ratio)
#
#     return(add_diff_var_ratio)
#   }
# )

# diff_var_ratio_vec output:
# 0.9965733 1.0468751 1.1011141 1.1596511 1.2229863 1.2916743 1.3662580 1.4474317
# 1.5359567 1.6326739 1.7386619
