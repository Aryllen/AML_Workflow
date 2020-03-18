###################################################################
#' @title Multivariate and Univariate Modeling
#'
#' @author Nicole Kauer
#'
#' @description
#' Does multivariate and univariate modeling for comparison.
#'
#' Execute from terminal:
#'     `Rscript multi_uni_analysis.R <directory>`
#' 
#' Inputs:
#'   `directory`: directory with PatientsVsDrugs.csv and
#'     PatientsVsVariants.csv
#' Outputs:
#'   
#'   
#' Packages:
#'   Biobase, leaps, BMA, iterativeBMA, fs, purrr, HDCI, glmnet,
#'   networkBMA, tibble, RColorBrewer, ComplexHeatmap, ggplot2
###################################################################

library(Biobase)
library(leaps)
library(BMA)
library(iterativeBMA)
library(fs)
library(purrr)
library(HDCI)
library(glmnet)
library(networkBMA)
library(tibble)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)

# Get the directory folder for drug and variant data
# Note: Not paying accounting for input errors.
# Assuming input is exactly correct.
dir <- commandArgs(trailingOnly=TRUE)
dir <- fs::path(dir)

##-----------------------------------------##
## FUNCTIONS                               ##
##-----------------------------------------##

#' @title Get indices to train and test on
#'
#' @description Gets `n` indices via random sampling from non-`NA` labels
#' to use as the training set,
#' where `n` = floor( (# non-`NA` samples) * percent_train ).
#' Returns the label indices to include in the training and
#' testing sets.
#' 
#' @inheritParams ibma_train_test
#' @param percent_train The percentage of patients to train on;
#'   floor used if number of patients does not divide evenly;
#'   numeric between 0 and 1
get_train_test_indices <- function(drugs, percent_train, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  ## Need only non-NA labels
  non_na_indices <- which(!is.na(drugs))
  ## Number of labels to choose for training
  num_to_train <- floor(length(non_na_indices) * percent_train)
  ## Sample indices within range of non_na_indices
  to_train_indices <- sample(1:length(non_na_indices), num_to_train)
  ## Get indices for original label set
  to_train <- non_na_indices[to_train_indices]
  to_test <- non_na_indices[-to_train_indices]
  ## Send back info
  list(
    to_train = to_train,
    to_test = to_test
  )
}

#' @title Custom iterativeBMA train, test, pred
#'
#' @description Divides up data and labels into training and testing sets,
#' removes data columns (ie genes) that have all 0's in
#' training set. Trains model, and tests on test data. Returns
#' results, including multiple error scores and information
#' about the model and inputs to the model.
#'
#' @param variants Full variant matrix with patients as rows and
#'   genes as columns
#' @param drugs Binary vector labels for drug sensitivity
#' @param drug_name Name of drug
#' @param to_train Vector of indices indicating the patients/rows to
#'   train on
#' @param to_test Vector of indices indicating the patients/rows to
#'   train on
#' @return List with:
#'   model: the trained model
#'   pred: the test prediction
#'   brier_score: brier_score / # test labels
#'   brier_score_errors: brier_score with 0/1 predictions / # test labels
#'   baseline_error: 1 - "true" labels / # labels
#'   inputs: some input information:
#'     drug_name: 
#'     trained_on: indices that model was trained on
#'     tested_on: indices that model was tested on
ibma_train_test <- function(variants, drugs, drug_name, to_train, to_test) {
  train_data <- variants[to_train, ]
  ## Get labels for same patients in training data
  train_labels <- drugs[to_train]
  
  ## Train model
  training_results <- iterateBMAlm(
    train_data,
    train_labels
  )
  
  train_data_subset <- train_data[, training_results$namesx[training_results$probne0 > 0]]
  training_results_new <- bicreg(train_data_subset, train_labels)
  
  ## Get test data and labels
  test_data <- variants[to_test, training_results$namesx]
  test_labels <- drugs[to_test]
  
  ## Get predictions for test data
  pred <- predict(
    training_results_new,
    test_data
  )
  ## Send back info
  list(
    model = training_results,
    new_model = training_results_new,
    pred = pred,
    test_labels = test_labels,
    inputs = list(
      drug_name = drug_name,
      trained_on = rownames(variants)[to_train],
      tested_on = rownames(variants)[to_test]
    )
  )
}

#' @title Get probne0 table of Genes versus Drugs
#'
#' @description Gathers all the probne0 values for each gene-drug pair
#' in the final bicreg models.
#'
#' @param models List of model objects given by the custom
#'   ibma_train_test()
#' @param genes Full set of possible genes used
#' @param threshold The minimum value, between 0 and 100,
#'   of probne0 to be included
get_probne0_all_drugs <- function(models, genes, threshold) {
  all_probne0 <- tibble::tibble(Genes = genes)
  for (model in models) {
    if (!inherits(model, "try-error")) {
      probne0 <- tibble::tibble(
        genes = model$new_model$namesx[model$new_model$probne0 > threshold],
        value = model$new_model$probne0[model$new_model$probne0 > threshold]
      )
      colnames(probne0) <- c("Genes", model$inputs$drug_name)
      all_probne0 <- dplyr::left_join(all_probne0, probne0, by = "Genes")
    }
  }
  ## Remove all genes that were never chosen
  na_rows <- apply(all_probne0[, 2:dim(all_probne0)[2]], 1, function(x) {
    all(is.na(x))
  })
  if (any(na_rows)) {
    all_probne0 <- all_probne0[!na_rows, ] 
  }
  all_probne0
}

#' @title Do linear regression training and testing
#'
#' @description Divides data into training and testing sets. Creates linear
#' regression model with training data. Uses test data for prediction and
#' calculates pearson, spearman, and RSME. Also calculates RSME for a baseline
#' prediction of a slight jitter around the mean of the training data.
#' 
#' @param genes Vector of genes to include in model
#' @param drug_name Name of the drug
#' @param drugs Vector of drug EC50 values
#' @param variants Dataframe with all gene variants for all patients
#' @param to_train Vector of indices to train on
#' @param to_test vector of indices to test on
linreg_train_test <- function(genes, drug_name, drugs, variants, to_train, to_test) {
  train_var <- variants[to_train, genes, drop = FALSE]
  train_labels <- drugs[to_train]
  train_data <- cbind(train_var, EC50 = train_labels)
  formula <- glue::glue_collapse(genes, sep = " + ")
  formula <- glue::glue("EC50 ~ {formula}")
  lm_model <- lm(formula, data = train_data)
  
  test_var <- variants[to_test, genes, drop = FALSE]
  test_labels <- drugs[to_test]
  test_data <- cbind(test_var, EC50 = test_labels)
  pred <- predict.lm(lm_model, newdata = test_data)
  pearson <- cor(test_labels, pred, method = "pearson")
  spearman <- cor(test_labels, pred, method = "spearman")
  
  # RMSE
  RSS <- sum((pred - test_labels)^2)
  MSE <- RSS / length(test_labels)
  RSME <- sqrt(MSE)
  
  # Get a vector the length of the test labels that is a jittered mean of the
  # training labels
  fake_labels <- sapply(test_labels, function(x) {jitter(mean(train_labels))})
  RSS_fake <- sum((fake_labels - test_labels)^2)
  MSE_fake <- RSS_fake / length(test_labels)
  RSME_fake <- sqrt(MSE_fake)
  
  list(
    drug = drug_name,
    genes = genes,
    model = lm_model,
    inputs = list(trained_on = to_train, tested_on = to_test),
    prediction = pred,
    test_labels = test_labels,
    correlation = list(
      spearman = spearman,
      pearson = pearson,
      rsme = RSME,
      baseline_rsme = RSME_fake
    )
  )
}

#' @title Get model coefficients lm
#'
#' @desciption Gathers all the coefficients for each gene-drug pair
#' in the linear regression models.
#'
#' @param models List of model objects given by the custom
#'   ibma_train_test()
#' @param genes Full set of possible genes used
get_lm_coeff_all_drugs <- function(models, genes) {
  all_coeff <- tibble::tibble(Genes = c("Intercept)", genes))
  for (model in models) {
    if (!inherits(model, "try-error")) {
      coeff <- tibble::tibble(
        Genes = names(model$model$coefficients),
        value = model$model$coefficients
      )
      colnames(coeff) <- c("Genes", model$drug)
      all_coeff <- dplyr::left_join(all_coeff, coeff, by = "Genes")
    }
  }
  ## Remove all genes that were never chosen
  na_rows <- apply(all_coeff[, 2:dim(all_coeff)[2]], 1, function(x) {
    all(is.na(x))
  })
  if (any(na_rows)) {
    all_coeff <- all_coeff[!na_rows, ] 
  }
  all_coeff
}

#' @title Do both univariate and mulitvariate linear regression models
#'
#' @description Runs the linear regression modeling for univariate analysis and
#' multivariate analysis.
#'
#' @inheritParams linreg_train_test
#' @param uni_var univariate gene
#' @param multi_var multivariate gene vector
train_test_uni_multi <- function(variants, drugs, drug_name, to_train, to_test,
                                 uni_var, multi_var) {
  uni_model <- linreg_train_test(uni_var, drug_name, drugs, variants, to_train, to_test)
  multi_model_all <- linreg_train_test(multi_var, drug_name, drugs, variants, to_train, to_test)
  list(
    uni = uni_model,
    multi_all = multi_model_all
  )
}

##-----------------------------------------##

## Get data
var <- read.csv(
  fs::path(dir, "PatientsVsVariants.csv"),
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
drug <- read.csv(
  fs::path(dir, "PatientsVsDrugs.csv"),
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Fix variant names to not have any punctuation
new_var_names <- unlist(lapply(
  colnames(var),
  function(x) {
    gsub("[[:punct:]]+", "", x)
  }
))
colnames(var) <- new_var_names
## Make variants into 0/1 matrix
var_01 <- var
var_01[var > 0] <- 1
## Remove genes with < 5 or > 64 variants
keep_genes <- apply(var_01, 2, function(x){sum(x) > 4 && sum(x) < 65})
var_01_sub <- var_01[keep_genes]

# Fix drug names to not have any punctuation
new_drug_names <- unlist(lapply(
  colnames(drug),
  function(x) {
    gsub("[[:punct:]]+", "", x)
  }
))
colnames(drug) <- new_drug_names
# Only use drugs with data for all patients
no_na <- apply(drug, 2, function(x) {all(!is.na(x))})
drugs_to_try <- drug[, no_na]

# Do train/test/pred for all drugs for feature selection.
# This has the seed set to get the same results every time.
# Ideally, this would be run 5 times without a seed set and then the
# model results averaged.
model_results <- purrr::map2(drugs_to_try, colnames(drugs_to_try), function(x, y) {
  indices <- get_train_test_indices(x, 0.8, 42)
  try(
    ibma_train_test(var_01_sub, x, y, indices$to_train, indices$to_test) 
  )
})
probne0_table <- get_probne0_all_drugs(model_results, colnames(var_01_sub), 5)
write.csv(probne0_table, fs::path(dir, "probne0.csv"))

# Do linreg for uni & multi but for the specific drugs that had low p-values
drugs_oi <- c("PLX4720", "MLN8237", "GDC0449", "Vinblastine", "Pazopanib",
              "Gefitinib", "Vincristine", "Daunorubicin", "Masitinib",
              "Irinotecan", "E7080", "AT7519", "Bexarotene", "SNS032",
              "BIBF1120", "Thioguanine", "E7080", "PD0332991", "Masitinib",
              "PD0325901", "Pp242", "SNS032", "Cabazitaxel", "Olaparib",
              "MLN8237", "Etoposide", "Tipifarnib", "Temsirolimus", "BIBF1120",
              "Acrichine", "Acrichine")
uni_genes <- c("CBX5", "SUZ12", "STAG2", "CBX5", "CBX5", "SUZ12", "CBX5",
               "SUPT5H", "PTPN11", "PTPN11", "PTPN11", "NPM1", "GAS6", "NPM1",
               "PTPN11", "CBX5", "CBX5", "CBX5", "SUZ12", "NRAS", "NRAS",
               "TET1", "ABL1", "PTPN11", "CBL", "SUPT5H", "CBX5", "CBX5",
               "CBX5", "SF3A1", "CEP164")

# This is a helper function for getting the data out of the model results
results_single_drug <- function(models) {
  # Uni
  uni_var <- tibble::tribble(
    ~drug, ~Genes_uni, ~Spearman_uni, ~Pearson_uni, ~RSME_uni, ~baseline_RSME_uni,
    models$uni$drug, glue::glue_collapse(models$uni$genes, sep = ", "), models$uni$correlation$spearman, models$uni$correlation$pearson, models$uni$correlation$rsme, models$uni$correlation$baseline_rsme
  )
  if (length(models$multi_50) > 1) {
    multi_50 <- tibble::tribble(
      ~Genes_multi50, ~Spearman_multi50, ~Pearson_multi50, ~RSME_multi50, ~baseline_RSME_multi50,
      glue::glue_collapse(models$multi_50$genes, sep = ", "), models$multi_50$correlation$spearman, models$multi_50$correlation$pearson, models$multi_50$correlation$rsme, models$multi_50$correlation$baseline_rsme
    ) 
  } else {
    multi_50 <- tibble::tribble(
      ~Genes_multi50, ~Spearman_multi50, ~Pearson_multi50, ~RSME_multi50, ~baseline_RSME_multi50,
      NA, NA, NA, NA, NA
    )
  }
  cbind(uni_var, multi_50)
}

# This function is gross and needs to be replaced with a dynamic way of
# doing the univariate/multivariate modeling. The problem stems mostly from
# how the data is returned from the BMA for feature selection function.
uni_multi_all <- function() {
  # train/test on same sets
  indices <- get_train_test_indices(drugs_to_try[, drugs_oi[1]], 0.8)
  
  # This needs to be done dynamically, not statically like this.
  # Will need to change how the initial multivarite model is stored. The list
  # makes it difficult to reference the specific drug correctly.
  # plx4720
  multi_var_50 <- model_results$PLX4720$new_model$namesx[model_results$PLX4720$new_model$probne0 >= 50]
  plx4720 <- list(
    uni = try(linreg_train_test("CBX5", "PLX4720", drugs_to_try[, "PLX4720"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "PLX4720", drugs_to_try[, "PLX4720"], var_01_sub, indices$to_train, indices$to_test))
  )
  # MLN8237
  multi_var_50 <- model_results$MLN8237$new_model$namesx[model_results$MLN8237$new_model$probne0 >= 50]
  MLN8237_1 <- list(
    uni = try(linreg_train_test("SUZ12", "MLN8237", drugs_to_try[, "MLN8237"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "MLN8237", drugs_to_try[, "MLN8237"], var_01_sub, indices$to_train, indices$to_test))
  )
  # GDC0449
  multi_var_50 <- model_results$GDC0449$new_model$namesx[model_results$GDC0449$new_model$probne0 >= 50]
  GDC0449 <- list(
    uni = try(linreg_train_test("STAG2", "GDC0449", drugs_to_try[, "GDC0449"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "GDC0449", drugs_to_try[, "GDC0449"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Vinblastine
  multi_var_50 <- model_results$Vinblastine$new_model$namesx[model_results$Vinblastine$new_model$probne0 >= 50]
  Vinblastine <- list(
    uni = try(linreg_train_test("CBX5", "Vinblastine", drugs_to_try[, "Vinblastine"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Vinblastine", drugs_to_try[, "Vinblastine"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Pazopanib
  multi_var_50 <- model_results$Pazopanib$new_model$namesx[model_results$Pazopanib$new_model$probne0 >= 50]
  Pazopanib <- list(
    uni = try(linreg_train_test("CBX5", "Pazopanib", drugs_to_try[, "Pazopanib"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Pazopanib", drugs_to_try[, "Pazopanib"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Gefitinib
  multi_var_50 <- model_results$Gefitinib$new_model$namesx[model_results$Gefitinib$new_model$probne0 >= 50]
  Gefitinib <- list(
    uni = try(linreg_train_test("SUZ12", "Gefitinib", drugs_to_try[, "Gefitinib"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Gefitinib", drugs_to_try[, "Gefitinib"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Vincristine
  multi_var_50 <- model_results$Vincristine$new_model$namesx[model_results$Vincristine$new_model$probne0 >= 50]
  Vincristine <- list(
    uni = try(linreg_train_test("CBX5", "Vincristine", drugs_to_try[, "Vincristine"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Vincristine", drugs_to_try[, "Vincristine"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Daunorubicin
  multi_var_50 <- model_results$Daunorubicin$new_model$namesx[model_results$Daunorubicin$new_model$probne0 >= 50]
  Daunorubicin <- list(
    uni = try(linreg_train_test("SUPT5H", "Daunorubicin", drugs_to_try[, "Daunorubicin"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Daunorubicin", drugs_to_try[, "Daunorubicin"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Masitinib
  multi_var_50 <- model_results$Masitinib$new_model$namesx[model_results$Masitinib$new_model$probne0 >= 50]
  Masitinib_1 <- list(
    uni = try(linreg_train_test("PTPN11", "Masitinib", drugs_to_try[, "Masitinib"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Masitinib", drugs_to_try[, "Masitinib"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Irinotecan
  multi_var_50 <- model_results$Irinotecan$new_model$namesx[model_results$Irinotecan$new_model$probne0 >= 50]
  Irinotecan <- list(
    uni = try(linreg_train_test("PTPN11", "Irinotecan", drugs_to_try[, "Irinotecan"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Irinotecan", drugs_to_try[, "Irinotecan"], var_01_sub, indices$to_train, indices$to_test))
  )
  # E7080
  multi_var_50 <- model_results$E7080$new_model$namesx[model_results$E7080$new_model$probne0 >= 50]
  E7080_1 <- list(
    uni = try(linreg_train_test("PTPN11", "E7080", drugs_to_try[, "E7080"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "E7080", drugs_to_try[, "E7080"], var_01_sub, indices$to_train, indices$to_test))
  )
  # AT7519
  multi_var_50 <- model_results$AT7519$new_model$namesx[model_results$AT7519$new_model$probne0 >= 50]
  AT7519 <- list(
    uni = try(linreg_train_test("NPM1", "AT7519", drugs_to_try[, "AT7519"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "AT7519", drugs_to_try[, "AT7519"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Bexarotene
  multi_var_50 <- model_results$Bexarotene$new_model$namesx[model_results$Bexarotene$new_model$probne0 >= 50]
  Bexarotene <- list(
    uni = try(linreg_train_test("GAS6", "Bexarotene", drugs_to_try[, "Bexarotene"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Bexarotene", drugs_to_try[, "Bexarotene"], var_01_sub, indices$to_train, indices$to_test))
  )
  # SNS032
  multi_var_50 <- model_results$SNS032$new_model$namesx[model_results$SNS032$new_model$probne0 >= 50]
  SNS032_1 <- list(
    uni = try(linreg_train_test("NPM1", "SNS032", drugs_to_try[, "SNS032"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "SNS032", drugs_to_try[, "SNS032"], var_01_sub, indices$to_train, indices$to_test))
  )
  # BIBF1120
  multi_var_50 <- model_results$BIBF1120$new_model$namesx[model_results$BIBF1120$new_model$probne0 >= 50]
  BIBF1120_1 <- list(
    uni = try(linreg_train_test("PTPN11", "BIBF1120", drugs_to_try[, "BIBF1120"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "BIBF1120", drugs_to_try[, "BIBF1120"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Thioguanine
  multi_var_50 <- model_results$Thioguanine$new_model$namesx[model_results$Thioguanine$new_model$probne0 >= 50]
  Thioguanine <- list(
    uni = try(linreg_train_test("CBX5", "Thioguanine", drugs_to_try[, "Thioguanine"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Thioguanine", drugs_to_try[, "Thioguanine"], var_01_sub, indices$to_train, indices$to_test))
  )
  # E7080
  multi_var_50 <- model_results$E7080$new_model$namesx[model_results$E7080$new_model$probne0 >= 50]
  E7080_2 <- list(
    uni = try(linreg_train_test("CBX5", "E7080", drugs_to_try[, "E7080"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "E7080", drugs_to_try[, "E7080"], var_01_sub, indices$to_train, indices$to_test))
  )
  # PD0332991
  multi_var_50 <- model_results$PD0332991$new_model$namesx[model_results$PD0332991$new_model$probne0 >= 50]
  PD0332991 <- list(
    uni = try(linreg_train_test("CBX5", "PD0332991", drugs_to_try[, "PD0332991"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "PD0332991", drugs_to_try[, "PD0332991"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Masitinib
  multi_var_50 <- model_results$Masitinib$new_model$namesx[model_results$Masitinib$new_model$probne0 >= 50]
  Masitinib_2 <- list(
    uni = try(linreg_train_test("SUZ12", "Masitinib", drugs_to_try[, "Masitinib"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Masitinib", drugs_to_try[, "Masitinib"], var_01_sub, indices$to_train, indices$to_test))
  )
  # PD0325901
  multi_var_50 <- model_results$PD0325901$new_model$namesx[model_results$PD0325901$new_model$probne0 >= 50]
  PD0325901 <- list(
    uni = try(linreg_train_test("NRAS", "PD0325901", drugs_to_try[, "PD0325901"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "PD0325901", drugs_to_try[, "PD0325901"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Pp242
  multi_var_50 <- model_results$Pp242$new_model$namesx[model_results$Pp242$new_model$probne0 >= 50]
  Pp242 <- list(
    uni = try(linreg_train_test("NRAS", "Pp242", drugs_to_try[, "Pp242"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Pp242", drugs_to_try[, "Pp242"], var_01_sub, indices$to_train, indices$to_test))
  )
  # SNS032
  multi_var_50 <- model_results$SNS032$new_model$namesx[model_results$SNS032$new_model$probne0 >= 50]
  SNS032_2 <- list(
    uni = try(linreg_train_test("TET1", "SNS032", drugs_to_try[, "SNS032"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "SNS032", drugs_to_try[, "SNS032"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Olaparib
  multi_var_50 <- model_results$Olaparib$new_model$namesx[model_results$Olaparib$new_model$probne0 >= 50]
  Olaparib <- list(
    uni = try(linreg_train_test("PTPN11", "Olaparib", drugs_to_try[, "Olaparib"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Olaparib", drugs_to_try[, "Olaparib"], var_01_sub, indices$to_train, indices$to_test))
  )
  # MLN8237
  multi_var_50 <- model_results$MLN8237$new_model$namesx[model_results$MLN8237$new_model$probne0 >= 50]
  MLN8237_2 <- list(
    uni = try(linreg_train_test("CBL", "MLN8237", drugs_to_try[, "MLN8237"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "MLN8237", drugs_to_try[, "MLN8237"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Etoposide
  multi_var_50 <- model_results$Etoposide$new_model$namesx[model_results$Etoposide$new_model$probne0 >= 50]
  Etoposide <- list(
    uni = try(linreg_train_test("SUPT5H", "Etoposide", drugs_to_try[, "Etoposide"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Etoposide", drugs_to_try[, "Etoposide"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Tipifarnib
  multi_var_50 <- model_results$Tipifarnib$new_model$namesx[model_results$Tipifarnib$new_model$probne0 >= 50]
  Tipifarnib <- list(
    uni = try(linreg_train_test("CBX5", "Tipifarnib", drugs_to_try[, "Tipifarnib"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Tipifarnib", drugs_to_try[, "Tipifarnib"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Temsirolimus
  multi_var_50 <- model_results$Temsirolimus$new_model$namesx[model_results$Temsirolimus$new_model$probne0 >= 50]
  Temsirolimus <- list(
    uni = try(linreg_train_test("CBX5", "Temsirolimus", drugs_to_try[, "Temsirolimus"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Temsirolimus", drugs_to_try[, "Temsirolimus"], var_01_sub, indices$to_train, indices$to_test))
  )
  # BIBF1120
  multi_var_50 <- model_results$BIBF1120$new_model$namesx[model_results$BIBF1120$new_model$probne0 >= 50]
  BIBF1120_2 <- list(
    uni = try(linreg_train_test("CBX5", "BIBF1120", drugs_to_try[, "BIBF1120"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "BIBF1120", drugs_to_try[, "BIBF1120"], var_01_sub, indices$to_train, indices$to_test))
  )
  # Acrichine
  multi_var_50 <- model_results$Acrichine$new_model$namesx[model_results$Acrichine$new_model$probne0 >= 50]
  Acrichine_1 <- list(
    uni = try(linreg_train_test("SF3A1", "Acrichine", drugs_to_try[, "Acrichine"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Acrichine", drugs_to_try[, "Acrichine"], var_01_sub, indices$to_train, indices$to_test))
  )
  Acrichine_2 <- list(
    uni = try(linreg_train_test("CEP164", "Acrichine", drugs_to_try[, "Acrichine"], var_01_sub, indices$to_train, indices$to_test)),
    multi_50 = try(linreg_train_test(multi_var_50, "Acrichine", drugs_to_try[, "Acrichine"], var_01_sub, indices$to_train, indices$to_test))
  )
  
  # Results have _1 and _2 after drug names that are the same due to not being
  # able to have duplicate row/column names in R.
  all_results <- rbind(
    results_single_drug(plx4720),
    results_single_drug(MLN8237_1),
    results_single_drug(GDC0449),
    results_single_drug(Vinblastine),
    results_single_drug(Pazopanib),
    results_single_drug(Gefitinib),
    results_single_drug(Vincristine),
    results_single_drug(Daunorubicin),
    results_single_drug(Masitinib_1),
    results_single_drug(Irinotecan),
    results_single_drug(E7080_1),
    results_single_drug(AT7519),
    results_single_drug(Bexarotene),
    results_single_drug(SNS032_1),
    results_single_drug(BIBF1120_1),
    results_single_drug(Thioguanine),
    results_single_drug(E7080_2),
    results_single_drug(PD0332991),
    results_single_drug(Masitinib_2),
    results_single_drug(PD0325901),
    results_single_drug(Pp242),
    results_single_drug(SNS032_2),
    results_single_drug(Olaparib),
    results_single_drug(MLN8237_2),
    results_single_drug(Etoposide),
    results_single_drug(Tipifarnib),
    results_single_drug(Temsirolimus),
    results_single_drug(BIBF1120_2),
    results_single_drug(Acrichine_1),
    results_single_drug(Acrichine_2)
  )
}

# Do 5 times and average. This is another area for code improvement...
# Also, this won't include the standard deviation when I do the average
res_1 <- uni_multi_all()
res_2 <- uni_multi_all()
res_3 <- uni_multi_all()
res_4 <- uni_multi_all()
res_5 <- uni_multi_all()

# Need to make NA == 0 or will not sum correctly
res_1[is.na(res_1)] <- 0
res_2[is.na(res_2)] <- 0
res_3[is.na(res_3)] <- 0
res_4[is.na(res_4)] <- 0
res_5[is.na(res_5)] <- 0

numeric_cols <- c(3, 4, 5, 6, 8, 9, 10, 11)
avg_res_numeric <- (res_1[, numeric_cols] + res_2[, numeric_cols] + res_3[, numeric_cols] + res_4[, numeric_cols] + res_5[, numeric_cols])/5
avg_res <- res_1 # Note that the multivariate genes and univariate genes do not change; just numeric values
avg_res[, c(3, 4, 5, 6)] <- avg_res_numeric[, c(1, 2, 3, 4)]
avg_res[, c(8, 9, 10, 11)] <- avg_res_numeric[, c(5, 6, 7, 8)]
write.csv(avg_res, fs::path(dir, "average_model_results.csv"))
