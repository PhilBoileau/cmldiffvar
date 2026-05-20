test_that(paste(
  "melt_tte_data_fun() melts the time-to-event data with correct number of",
  "repeated measures"
  ),
  {
    set.seed(18341)

    # generate the data
    sample_size <- 15
    wide_tte_dat_tbl <- generate_test_data(n_obs = 15)

    # set the restriction time
    time_cutoff <- 100

    # melt the data
    long_tte_dat_tbl <- melt_tte_data_fun(
      wide_data_tbl = wide_tte_dat_tbl,
      baseline_var_names = "w",
      treatment_var_name = "a",
      outcome_var_name = "time",
      censoring_var_name = "censoring",
      time_cutoff = time_cutoff
    )

    # confirm that there are 15 unique observations in the long data
    num_unique_ids <- long_tte_dat_tbl |>
      pull(cmldiffvar_id) |>
      unique() |>
      length()
    expect_equal(sample_size, num_unique_ids)

    # confirm that each observation has one row for each of the possible time
    # entries up to the cutoff
    expect_equal(time_cutoff * sample_size, nrow(long_tte_dat_tbl))

    # confirm that the nuisance indicator only contains rows associated with
    # any event times prior to the time_cutoff
    unique_times <- wide_tte_dat_tbl |>
      filter(time <= time_cutoff) |>
      pull(time) |>
      unique() |>
      length()
    nuisance_estimation_tbl <- long_tte_dat_tbl |>
      filter(cmldiffvar_for_nuisance_estimation == 1)
    expect_equal(unique_times * sample_size, nrow(nuisance_estimation_tbl))

  }
)
