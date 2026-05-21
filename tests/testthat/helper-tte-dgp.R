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
    censoring = censoring
  )

  return(wide_tbl)
}
