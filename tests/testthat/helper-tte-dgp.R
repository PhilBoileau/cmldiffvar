library(dplyr)

# function for randomly generating a basic TTE outcome dataset for testing
generate_test_data <- function(
    n_obs = 200,
    null_marginal = FALSE,
    hazard_model = TRUE
) {

  # baseline covariate
  w <- rbinom(n = n_obs, size = 1, prob = 0.5)

  # exposure
  prop_score <- plogis(-0.2 + 0.7 * w)
  a <- rbinom(n_obs, 1, prop_score)

  # generate data according to hazards model
  if (hazard_model) {
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
      cond_surv_hazard <- function(time, exposure, w) {
        (time < max_time) / (75 + exp(2 * (exposure + w + w * exposure))) +
          (time == max_time)
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
  } else {

    # otherwise use a simple AFT model
    eps_0 <- rnorm(n_obs)
    censor_time_0 <- floor(rexp(n_obs, rate = 0.01))
    censor_time_1 <- floor(rexp(n_obs, rate = 0.01))
    if (null_marginal) {
      failure_time_1 <- failure_time_0 <- floor(exp(3.5 - 0.4 * w + eps_0))
    } else {
      failure_time_0 <- floor(exp(3.5 - 0.4 * w + eps_0))
      failure_time_1 <- failure_time_0 * 1.1
    }

  }

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

# Hazard DGP ----

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
# add_diff_var_ratio_tbl <- lapply(
#   seq(from = 0, to = 5),
#   function(trunc_value) {
#
#     # compute the additive diff var ratio under the homogeneous additive
#     # treatment effect assumption
#     add_diff_var_ratio_tbl <- pop_sim_tbl |>
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
#       mutate(trunc_value = trunc_value)
#
#     return(add_diff_var_ratio_tbl)
#   }
# ) |>
#   bind_rows()

# diff_var_ratio_tbl output:
#   var_0 var_1 add_diff_var_ratio trunc_value
#   <dbl> <dbl>              <dbl>       <int>
# 1  280.  279.              0.997           0
# 2  266.  279.              1.05            1
# 3  253.  279.              1.10            2
# 4  241.  279.              1.16            3
# 5  228.  279.              1.22            4
# 6  216.  279.              1.29            5


# AFT DGP ----

## Null marginal model ----

# # approximate relevant parameter values at time_cutoff = 50
# set.seed(486415)
# pop_sim_tbl <- generate_test_data(
#   n_obs = 100000, null_marginal = TRUE, hazard_model = FALSE
# )
# time_cutoff <- 100
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
#     rmst_ratio = rmst_1 / rmst_0,
#     diff_var_ratio = var_1 / var_0
#   )
#
# # recorded values
# # $ rmst_0         <dbl> 47.24593
# # $ rmst_1         <dbl> 47.05074
# # $ var_0          <dbl> 1026.412
# # $ var_1          <dbl> 1030.285
# # $ rmst_ratio     <dbl> 0.9958686
# # $ diff_var_ratio <dbl> 1.003774
#
# # compute the ratios of truncated survival times
# mult_diff_var_ratio_tbl <- lapply(
#   seq(from = 0.90, to = 1.0, by = 0.01),
#   function(trunc_value) {
#
#     # compute the multiplicative diff var ratio under the homogeneous
#     # multiplicative treatment effect assumption
#     mult_diff_var_ratio_tbl <- pop_sim_tbl |>
#       mutate(
#         trunc_trunc_potential_time_0 = if_else(
#           potential_time_0 > (time_cutoff * trunc_value),
#           (time_cutoff * trunc_value), potential_time_0
#         ),
#         trunc_potential_time_1 = if_else(
#           potential_time_1 > time_cutoff, time_cutoff, potential_time_1
#         )
#       ) |>
#       summarise(
#         var_0 = var(trunc_trunc_potential_time_0),
#         var_1 = var(trunc_potential_time_1),
#         mult_diff_var_ratio = var_1 / var_0
#       ) |>
#       mutate(trunc_value = trunc_value)
#
#     return(mult_diff_var_ratio_tbl)
#   }
# ) |>
#   bind_rows()
#
# # mult_diff_var_ratio_tbl output
# #    var_0 var_1 mult_diff_var_ratio trunc_value
# #    <dbl> <dbl>               <dbl>       <dbl>
# #  1  841. 1030.                1.23        0.9
# #  2  859. 1030.                1.20        0.91
# #  3  877. 1030.                1.17        0.92
# #  4  896. 1030.                1.15        0.93
# #  5  914. 1030.                1.13        0.94
# #  6  933. 1030.                1.10        0.95
# #  7  951. 1030.                1.08        0.96
# #  8  970. 1030.                1.06        0.97
# #  9  989. 1030.                1.04        0.98
# # 10 1007. 1030.                1.02        0.99
# # 11 1026. 1030.                1.00        1


## Heterogenous model ----

# approximate relevant parameter values at time_cutoff = 100
set.seed(19984)
pop_sim_tbl <- generate_test_data(
  n_obs = 100000, null_marginal = TRUE, hazard_model = FALSE
)
time_cutoff <- 100
pop_params_tbl <- pop_sim_tbl |>
  mutate(
    trunc_potential_time_0 = if_else(
      potential_time_0 > time_cutoff, time_cutoff, potential_time_0
    ),
    trunc_potential_time_1 = if_else(
      potential_time_1 > time_cutoff, time_cutoff, potential_time_1
    )
  ) |>
  summarise(
    rmst_0 = mean(trunc_potential_time_0),
    rmst_1 = mean(trunc_potential_time_1),
    var_0 = var(trunc_potential_time_0),
    var_1 = var(trunc_potential_time_1),
    rmst_ratio = rmst_1 / rmst_0,
    diff_var_ratio = var_1 / var_0
  )

# recorded values
# $ rmst_0         <dbl> 36.92097
# $ rmst_1         <dbl> 36.92097
# $ var_0          <dbl> 916.7922
# $ var_1          <dbl> 916.7922
# $ rmst_ratio     <dbl> 1
# $ diff_var_ratio <dbl> 1

# compute the ratios of truncated survival times
mult_diff_var_ratio_tbl <- lapply(
  seq(from = 0.90, to = 1.00, by = 0.01),
  function(trunc_value) {

    # compute the multiplicative diff var ratio under the homogeneous
    # multiplicative treatment effect assumption
    mult_diff_var_ratio_tbl <- pop_sim_tbl |>
      mutate(
        trunc_trunc_potential_time_0 = if_else(
          potential_time_0 > (time_cutoff * trunc_value),
          (time_cutoff * trunc_value), potential_time_0
        ),
        trunc_potential_time_1 = if_else(
          potential_time_1 > time_cutoff, time_cutoff, potential_time_1
        )
      ) |>
      summarise(
        var_0 = var(trunc_trunc_potential_time_0),
        var_1 = var(trunc_potential_time_1),
        mult_diff_var_ratio = var_1 / var_0,
        standardized_diff_var = mult_diff_var_ratio - (1 / trunc_value)^2
      ) |>
      mutate(trunc_value = trunc_value)

    return(mult_diff_var_ratio_tbl)
  }
) |>
  bind_rows()

# mult_diff_var_ratio_tbl output
#    var_0 var_1 mult_diff_var_ratio standardized_diff_var trunc_value
#    <dbl> <dbl>               <dbl>                 <dbl>       <dbl>
#  1  790.  917.                1.16              -0.0742         0.9
#  2  803.  917.                1.14              -0.0656         0.91
#  3  816.  917.                1.12              -0.0574         0.92
#  4  828.  917.                1.11              -0.0494         0.93
#  5  841.  917.                1.09              -0.0417         0.94
#  6  854.  917.                1.07              -0.0342         0.95
#  7  866.  917.                1.06              -0.0270         0.96
#  8  879.  917.                1.04              -0.0199         0.97
#  9  892.  917.                1.03              -0.0131         0.98
# 10  904.  917.                1.01              -0.00646        0.99
# 11  917.  917.                1                  0              1
