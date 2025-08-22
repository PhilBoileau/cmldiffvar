test_that("serial cmldiffvarr estimates diff vars without errors", {

  library(dplyr)
  library(SuperLearner)

  set.seed(2352314)

  # sample tibble
  sample_tbl <- slice_sample(toy_population_tbl, n = 500)

  # TMLE for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
    )

  # TMLE for absolute diff var with known propensity scores
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_var_name = "propensity_score",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # TMLE for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # cross-fitted TMLE for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # cross-fitted TMLE for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # cross-fitted TMLE for relative diff var with known propensity scores
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_var_name = "propensity_score",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # one-step for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # one-step for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # one-step for relative diff var with known propensity score
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_var_name = "propensity_score",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # cross-fitted one-step for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # cross-fitted one-step for absolute diff var with known propensity scores
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        propensity_score_var_name = "propensity_score",
        outcome_var_name = "outcome",
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

  # cross-fitted one-step for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )

})

test_that("parallelized cmldiffvarr estimates diff vars without errors",{

  library(dplyr)
  library(future)
  library(SuperLearner)

  plan(sequential)

  set.seed(2352314)

  # sample tibble
  sample_tbl <- slice_sample(toy_population_tbl, n = 500)

  # TMLE for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )

  # TMLE for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )

  # cross-fitted TMLE for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )

  # cross-fitted TMLE for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )

  # one-step for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )

  # one-step for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )

  # cross-fitted one-step for absolute diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )

  # cross-fitted one-step for relative diff var
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "relative",
        estimator_type = "one-step",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.mean", "SL.glm.gamma.identity", "SL.gam.gamma.log"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = TRUE,
        num_cross_fit_folds = 5,
        parallel = TRUE
      )
  )
})

test_that("SuperLearner correctly calls predict.SL.torch
          when using torch-based learners",{

  library(dplyr)
  library(SuperLearner)

  set.seed(2352314)

  # sample tibble
  sample_tbl <- slice_sample(toy_population_tbl, n = 500)

  # cond_exp_sq_outcome is predicted using SL.nnet.torch.softplus
  expect_no_error(
    sample_tbl |>
      cmldiffvar(
        estimand_type = "absolute",
        estimator_type = "tmle",
        confidence_level = 0.95,
        confounder_var_names = "confounder",
        treatment_var_name = "treatment",
        outcome_var_name = "outcome",
        propensity_score_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_outcome_library = c("SL.mean", "SL.glm", "SL.earth"),
        cond_exp_sq_outcome_library = c(
          "SL.nnet.torch.softplus"
        ),
        num_nuisance_sl_folds = 5,
        cross_fit = FALSE,
        num_cross_fit_folds = 5,
        parallel = FALSE
      )
  )
})

