###################################################################
#' @title Fix drug synonyms
#'
#' @author Nicole Kauer
#'
#' @description
#' Removes synonyms from the Compounds_Remove.csv. Assumes an
#' output directory has been created in the drug sensitivity
#' directory given, and that a Compound_Synonyms.csv file exists
#' in that directory.
#'
#' Execute from terminal:
#'     `Rscript DrugSensitivity_FixSynonymDrugs.R <directory>`
#' 
#' Outputs:
#'   Compounds_Remove.csv: Overwrites original Compounds_Remove.csv
#'     with the synonym drugs removed.
###################################################################

# Line below is simply for RStudio. Turning off diagnostics because keeps giving
# annoying, useless warnings; apparently, this is a known issue.
# !diagnostics off

# For string matching
library("stringr", quietly=TRUE)
library("fs", quietly=TRUE)

# Get the directory folder for the drug sensitivity data
# Note: Not paying accounting for input errors. Assuming input is exactly correct.
drug.dir <- commandArgs(trailingOnly=TRUE)
drug.dir <- fs::path(drug.dir)

# The output directory is in the same input directory
output.dir <- fs::path(drug.dir, "Output")

# Open Compounds_Remove.csv for editing
compounds.remove <- read.csv(fs::path(output.dir, "Compounds_Remove.csv"), header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
compounds.all <- compounds.remove$Compound

# Get information on synonym file since can't open if empty
syn.info <- file.info(fs::path(output.dir, "Compound_Synonyms.csv"))
# Make sure synonym file is not empty; R reads this in different ways sometimes, hence the multiple conditions
if (!is.na(syn.info$size) & syn.info$size > 0) {
  # Open synonym file
  synonyms <- read.csv(fs::path(output.dir, "Compound_Synonyms.csv"), header = FALSE, na.strings = c("", " ", NA), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  
  # Add all compounds in first column of synonym list
  #    This allows for names that users want to be shorter but never appear in compound list in files
  compounds.all <- append(compounds.all, synonyms[,1])
  #apply(synonyms, 1, function (x) {compounds.all <- rbind(compounds.all, data.frame(Compound = x[1], Include = 1))})
  # Get the unique set in case synonym file added a double
  # NOTE: need to check that this works if one of the Include values is 0
  compounds.all <- unique(compounds.all, ignore.case = TRUE)
  # Transform compounds to be lowercase and have no punctuation
  # New vector of compounds so the original names stay normal
  #compounds.all.t <- gsub("[[:punct:] ]+", "", compounds.all)
  compounds.all.t <- tolower(compounds.all)
  # List of indices to remove
  comp.ind <- vector()
  
  # Remove all synonym compounds; assume first column in synonyms file is primary name
  # Go through each row, ignoring first column in the row and any NA values
  for (i in 1:nrow(synonyms)) {
    # Go through each synonym and check the location of the name in the full compounds set
    for (synonym in synonyms[i,-1]) {
      # Check if the synonym is NA
      if (!is.na(synonym)) {
        # Transform synonym the same way as the compounds
        #synonym.t <- gsub("[[:punct:] ]+", "", synonym)
        synonym.t <- tolower(synonym)
        # Get index for where name appears in full compound set
        syn.index <- grep(synonym.t, compounds.all.t, fixed = TRUE)
        # Check that an appropriate index was actually found (i.e. name exists in list)
        if (length(syn.index) > 0) {
          # Add to list for removal
          comp.ind <- append(comp.ind, syn.index)
        }  
      }
    }
  } 
} else {
  quit(save = "no")
}

if (length(comp.ind) > 0) {
  compounds.all <- compounds.all[-comp.ind] 
}

# Sort
compounds.all <- sort(compounds.all)

all.compounds.included <- data.frame(Compound = compounds.all, Include = 1)
# Set Include to be 0 where it was previously
# Get the genes that should be removed (i.e. Include = 0)
compounds.to.remove <- compounds.remove$Compound[compounds.remove$Include == 0]
if (length(compounds.remove)[1] > 0) {
  # Set genes to remove that were previously set
  for (drug in compounds.to.remove) {
    drug.ind <- which(all.compounds.included$Compound == drug)
    if (length(drug.ind) > 0) {
      all.compounds.included$Include[drug.ind] <- 0 
    }
  }
}
write.csv(all.compounds.included, file = fs::path(output.dir, "Compounds_Remove.csv"), row.names = FALSE)
