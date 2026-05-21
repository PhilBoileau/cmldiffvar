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
        (time < max_time) / (75 + exp(1 + 2 * w + 3 * exposure)) +
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
    eps_1 <- rnorm(n_obs)
    censor_time_0 <- floor(rexp(n_obs, rate = 0.01))
    censor_time_1 <- floor(rexp(n_obs, rate = 0.01))
    if (null_marginal) {
      failure_time_0 <- floor(exp(3.5 - 0.4 * w + eps_0))
      failure_time_1 <- floor(failure_time_0 * 1.1)
    } else {
      failure_time_0 <- floor(exp(3.5 - 0.4 * w + eps_0))
      failure_time_1 <- floor(exp(3.5 + 0.4 * w + 0.25 * a + eps_1))
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

## Homogeneous additive model with null treatment effect ----

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



## Heterogeneous additive model with null marginal RMST difference ----

# approximate relevant parameter values at time_cutoff = 50
# set.seed(948)
# pop_sim_tbl <- generate_test_data(n_obs = 10000, null_marginal = FALSE)
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
#   # group_by(w) |>
#   summarise(
#     rmst_0 = mean(trunc_potential_time_0),
#     rmst_1 = mean(trunc_potential_time_1),
#     var_0 = var(trunc_potential_time_0),
#     var_1 = var(trunc_potential_time_1),
#     rmst_diff = rmst_1 - rmst_0,
#     diff_var_ratio = var_1 / var_0
#   )

# recorded values
# $ rmst_0         <dbl> 59.0448
# $ rmst_1         <dbl> 80.0494
# $ var_0          <dbl> 1285.684
# $ var_1          <dbl> 995.0123
# $ rmst_diff      <dbl> 21.0046
# $ diff_var_ratio <dbl> 0.7739167

# compute the ratios of truncated survival times
# add_diff_var_ratio_tbl <- lapply(
#   seq(from = 15, to = 30),
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
#    <dbl> <dbl>              <dbl>       <int>
#  1  924.  995.               1.08          15
#  2  901.  995.               1.10          16
#  3  878.  995.               1.13          17
#  4  856.  995.               1.16          18
#  5  834.  995.               1.19          19
#  6  812.  995.               1.23          20
#  7  790.  995.               1.26          21
#  8  768.  995.               1.29          22
#  9  747.  995.               1.33          23
# 10  726.  995.               1.37          24
# 11  705.  995.               1.41          25
# 12  685.  995.               1.45          26
# 13  664.  995.               1.50          27
# 14  644.  995.               1.54          28
# 15  624.  995.               1.59          29
# 16  605.  995.               1.65          30

# AFT DGP ----

## Homogeneous non-null treatment effect model ----

# approximate relevant parameter values at time_cutoff = 100
# set.seed(19984)
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

# recorded values
# $ rmst_0         <dbl> 36.92097
# $ rmst_1         <dbl> 39.53811
# $ var_0          <dbl> 916.7922
# $ var_1          <dbl> 970.0174
# $ rmst_ratio     <dbl> 1.070885
# $ diff_var_ratio <dbl> 1.058056

# compute the ratios of truncated survival times
# mult_diff_var_ratio_tbl <- lapply(
#   seq(from = 0.85, to = 1.00, by = 0.01),
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
#         mult_diff_var_ratio = var_1 / var_0,
#         standardized_diff_var = mult_diff_var_ratio - (1 / trunc_value)^2
#       ) |>
#       mutate(trunc_value = trunc_value)
#
#     return(mult_diff_var_ratio_tbl)
#   }
# ) |>
#   bind_rows()

# mult_diff_var_ratio_tbl output
#    var_0 var_1 mult_diff_var_ratio standardized_diff_var trunc_value
#    <dbl> <dbl>               <dbl>                 <dbl>       <dbl>
#  1  726.  970.                1.34             -0.0481          0.85
#  2  739.  970.                1.31             -0.0393          0.86
#  3  752.  970.                1.29             -0.0308          0.87
#  4  765.  970.                1.27             -0.0226          0.88
#  5  777.  970.                1.25             -0.0146          0.89
#  6  790.  970.                1.23             -0.00681         0.9
#  7  803.  970.                1.21              0.000655        0.91
#  8  816.  970.                1.19              0.00787         0.92
#  9  828.  970.                1.17              0.0149          0.93
# 10  841.  970.                1.15              0.0216          0.94
# 11  854.  970.                1.14              0.0281          0.95
# 12  866.  970.                1.12              0.0344          0.96
# 13  879.  970.                1.10              0.0406          0.97
# 14  892.  970.                1.09              0.0466          0.98
# 15  904.  970.                1.07              0.0524          0.99
# 16  917.  970.                1.06              0.0581          1

# as expected, the standardized differential variance is approximately equal to
# 0 when the truncation value is near (1/1.1), indicating a homogeneous but
# non-null effect

## Heterogeneous treatment effect model ----
# approximate relevant parameter values at time_cutoff = 100
set.seed(824563)
pop_sim_tbl <- generate_test_data(
  n_obs = 100000, null_marginal = FALSE, hazard_model = FALSE
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
# $ rmst_0         <dbl> 37.12697
# $ rmst_1         <dbl> 52.82716
# $ var_0          <dbl> 917.5128
# $ var_1          <dbl> 1135.992
# $ rmst_ratio     <dbl> 1.422878
# $ diff_var_ratio <dbl> 1.238121

# compute the ratios of truncated survival times
mult_diff_var_ratio_tbl <- lapply(
  seq(from = 0.6, to = 0.75, by = 0.01),
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

# mult_diff_var_ratio_tbl output:
#    var_0 var_1 mult_diff_var_ratio standardized_diff_var trunc_value
#    <dbl> <dbl>               <dbl>                 <dbl>       <dbl>
#  1  406. 1136.                2.80                0.0222        0.6
#  2  418. 1136.                2.72                0.0290        0.61
#  3  431. 1136.                2.64                0.0361        0.62
#  4  443. 1136.                2.56                0.0434        0.63
#  5  456. 1136.                2.49                0.0506        0.64
#  6  469. 1136.                2.42                0.0578        0.65
#  7  481. 1136.                2.36                0.0651        0.66
#  8  494. 1136.                2.30                0.0723        0.67
#  9  507. 1136.                2.24                0.0794        0.68
# 10  519. 1136.                2.19                0.0865        0.69
# 11  532. 1136.                2.13                0.0935        0.7
# 12  545. 1136.                2.08                0.100         0.71
# 13  558. 1136.                2.04                0.107         0.72
# 14  571. 1136.                1.99                0.113         0.73
# 15  584. 1136.                1.95                0.120         0.74
# 16  597. 1136.                1.90                0.126         0.75
