###################################################################
#' @title Drug Sensitivity and Variant Gathering
#'
#' @author Nicole Kauer
#'
#' @description
#' Gathers variants and drug sensitivity data based on filters
#' and outputs csv files, one for variants and one for drug
#' sensitivity. Requires variant file directory, drug
#' sensitivity directory, and an output directory. Note that the
#' patients output in the final dataset are only those that
#' have both variant and drug sensitivity data.
#' 
#' For options, if no SIFT or polyphen filters chosen, then all
#' variants in dataset are included.
#' 
#' Execute from terminal:
#'   `Rscript Variant_Source_to_CSV.r -v <var_dir> -d <ds_dir> -o <output_dir>`
#'    Include any data extra options.
#'    
#' Options:
#'   `-v`, `--var.dir`: Directory for variant files
#'   `-d`, `--drug.dir`: Directory for drug sensitivity files
#'   `-o`, `--output.dir`: Directory to output files
#'   `-f`, `--frequency``): Minimum variant allele frequency to be included in dataset;
#'      default is 2.5%
#'   `-a`, `--alpha`): The alpha value for statistic analysis; default is 0.05
#'   `--sift.deleterious`: True to include variants with SIFT = deleterious;
#'      default is FALSE
#'   `--sift.tolerated`: True to include missense variants with SIFT = tolerated;
#'      default is FALSE
#'   `--poly.prob.damaging`: True to include missense variants with Polyphen = probably damaging;
#'      default is FALSE
#'   `--poly.poss.damaging`: True to nclude missense variants with Polyphen = possibly damaging;
#'      default is FALSE
#'   `--poly.benign`: True to include missense variants with Polyphen = benign;
#'      default is FALSE
#'   `--active.maximum`: Max EC50 value to be considered active; values greater are considered inactive;
#'      default is 1E-4 M
#'   `--active.minimum`: Min EC50 value to be considered highly active;
#'      default is 1E-12 M
#'
#' Outputs:
#'   PatientsVsVariants.csv: spreadsheet with counts of gene variants for
#'     all patients.
#'   PatientsVsDrugs.csv: spreadsheet with -log(EC50) values for all drugs
#'     and patients.
#'   Patient_subset: Two column spreadsheet. Column `Patient` has all
#'     patients in current dataset. Column `Include` has a 1 to keep
#'     patient in dataset; 0 otherwise.
#'   Compound_subset: Two column spreadsheet. Column `Compound` has all
#'     compounds in current dataset. Column `Include` has a 1 to keep
#'     compound in dataset; 0 otherwise.
###################################################################

library("optparse", quietly=TRUE)
library("readxl", quietly=TRUE)
library("stringr", quietly=TRUE)
library("fs", quietly=TRUE)

# Command line options
options <- list(
  make_option(c("-v", "--var.dir"), type="character", default=NULL, help="Directory for variant files."),
  make_option(c("-d", "--drug.dir"), type="character", default=NULL, help="Directory for drug sensitivity files."),
  make_option(c("-o", "--output.dir"), type="character", default=NULL, help="Directory to output files."),
  make_option(c("-f", "--frequency"), type="character", default="2.5", help="Minimum variant frequency to be included in dataset."),
  make_option(c("-a", "--alpha"), type="character", default="0.05", help="The alpha value for statistic analysis."),
  make_option(c("--sift.deleterious"), action = "store_true", default=FALSE, help="Include variants with SIFT = deleterious."),
  make_option(c("--sift.tolerated"), action = "store_true", default=FALSE, help="Include missense variants with SIFT = tolerated."),
  make_option(c("--poly.prob.damaging"), action = "store_true", default=FALSE, help="Include missense variants with Polyphen = probably damaging."),
  make_option(c("--poly.poss.damaging"), action = "store_true", default=FALSE, help="Include missense variants with Polyphen = possibly damaging."),
  make_option(c("--poly.benign"), action = "store_true", default=FALSE, help="Include missense variants with Polyphen = benign."),
  make_option(c("--active.maximum"), type="character", default="1E-4", help="Maximum EC50 value to be considered active; all values greater are considered inactive."),
  make_option(c("--active.minimum"), type="character", default="1E-12", help="Minimum EC50 value to be considered highly active; all values smaller are considered highly active."),
  make_option(c("-p", "--patient.subset"), action = "store_true", default=FALSE, help="Choose the patient subset, else the previously chosen subset is used."),
  make_option(c("-c", "--compound.subset"), action = "store_true", default=FALSE, help="Choose the drug compound subset, else the previously chosen subset is used.")
)

# Argument retriever
opt_parser <- OptionParser(option_list=options)
# Get arguments
opts <- parse_args(opt_parser)

# Make all the numeric values be numeric
opts$frequency <- as.numeric(opts$frequency)
opts$alpha <- as.numeric(opts$alpha)
opts$active.maximum <- as.numeric(opts$active.maximum)
opts$active.minimum <- as.numeric(opts$active.minimum)

# Need to have the directories input; verify that they were.
# Note: not checking that they are valid inputs; also, if I can just move the names along from the previous widget in
# BwB, then this should be unnecessary to check.
if (is.null(opts$var.dir)) {
  stop("No folder for the variant files was detected. Please input the directory in the command line option '-v <directory>'.")
} else {
  # Append /Output/ to end because the code needs it for the directory to map correctly to cleaned files
  opts$var.dir <- fs::path(opts$var.dir, "Output")
}
if (is.null(opts$drug.dir)) {
  stop("No folder for the drug sensitivity files was detected. Please input the directory in the command line option '-d <directory>'.")
} else {
  # Append / to end because the code needs it for the directory to map correctly to files cleaned files
  opts$drug.dir <- fs::path(opts$drug.dir, "Output")
}
if (is.null(opts$output.dir)) {
  stop("No folder for the output files was detected. Please input the directory in the command line option '-o <directory>'.")
}

# Output directories for variant and drug sensitivity data
clean.var.dir <- fs::path(opts$var.dir,"_Clean")
clean.drug.dir <- fs::path(opts$drug.dir,"_Clean")

# Need to make sure all directories exist before continuing
if (!fs::dir_exists(opts$output.dir)) {
  fs::dir_create(opts$output.dir)
}
if (!fs::dir_exists(clean.var.dir)) {
  stop("No clean variant data directory found")
}
if (!fs::dir_exists(clean.drug.dir)) {
  stop("No clean drug data directory found")
}

##### Drug Sensitivity Data Gathering #####

# Get lists of the files in the directories
filenames.drug.excel <- list.files(clean.drug.dir, pattern = "*.xlsx")
filenames.drug.csv <- list.files(clean.drug.dir, pattern = "*.csv")

# Patient list from drug sensitivity cleaning
patients.drug <- as.numeric(readLines(fs::path(opts$drug.dir, "Patients_DrugSensitivity.txt")))

# Compound list from Compounds_Remove.csv
compounds.remove <- read.csv(fs::path(opts$drug.dir, "Compounds_Remove.csv"), check.names = FALSE, header=TRUE, na.strings = c("", " ", NA), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
compounds <- compounds.remove$Compound[which(compounds.remove$Include == 1)]

# Synonyms from file
syn.info <- file.info(fs::path(opts$drug.dir, "Compound_Synonyms.csv"))
# Make sure synonym file is not empty; R reads this in different ways sometimes, hence the multiple conditions
if (!is.na(syn.info$size) & syn.info$size > 0) {
  synonyms <- read.csv(fs::path(opts$drug.dir, "Compound_Synonyms.csv"), check.names = FALSE, header = FALSE, na.strings = c("", " ", NA), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
} else {
  synonyms <- NULL
}

##### Function(s) #####
# Synonym row index function
# Input: drug.name (string) - name of drug
# Returns: row index for matching drug name in synonym file or 0 if not found
syns.row <- function (drug.name) {
  if (is.null(synonyms)) {
    return(0)
  }
  #drug <- gsub("[[:punct:] ]+", "", drug.name)
  drug <- tolower(drug.name)
  # Initialize to default return value if drug name not found
  match.index <- 0
  # Go through each row synonym file
  for (i in 1:nrow(synonyms)) {
    #synonyms.t <- gsub("[[:punct:] ]+", "", synonyms[i,])
    # Check if drug name is in the row; save index if found
    if (drug %in% (tolower(synonyms[i, ]))) {
      match.index <- i
    }
  }
  return(match.index)
}
#####################

# Dataframe for the collected EC50 drug sensitivity data
data.drug.ec50 <- data.frame(matrix(nrow=length(patients.drug), ncol=length(compounds)))
rownames(data.drug.ec50) <- patients.drug
colnames(data.drug.ec50) <- compounds
# Dataframe for the collected AUC drug sensitivity data
data.drug.auc <- data.frame(matrix(nrow=length(patients.drug), ncol=length(compounds)))
rownames(data.drug.auc) <- patients.drug
colnames(data.drug.auc) <- compounds

# Go through all Excel files and start collecting the drug sensitivity data
for (data.file in filenames.drug.excel) {

  # Get names of sheets in Excel workbook
  sheets <- excel_sheets(fs::path(clean.drug.dir,data.file))

  # Find the sheet that (hopefully) has static data based on previous sheets that (sometimes) has it
  # Assuming the first column with the name is the right one
  if (length(which(str_detect(sheets, coll("Fitted Parameters Static", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("Fitted Parameters Static", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("FitParameters", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("FitParameters", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("combined", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("combined", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("Fitted Parameters", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("Fitted Parameters", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("blast param statsort", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("blast param statsort", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(sheets, coll("Parameters static", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("Parameters static", ignore_case = TRUE)))[1]
  }  else if (length(which(str_detect(sheets, coll("4Database", ignore_case = TRUE)))) > 0) {
    sensitivity.sheet <- which(str_detect(sheets, coll("4Database", ignore_case = TRUE)))[1]
  }
  #### Should have a graceful end here and in other places where it could crash somehow...

  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  patient.index <- which(rownames(data.drug.ec50) == patient)

  # Open file; suppress unimportant messages, again
  drug.data <- suppressMessages(read_excel(fs::path(clean.drug.dir, data.file), sheet=sheets[sensitivity.sheet]))

  # Get column that has drug compounds
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  }

  # Get column that has either EC50 or IC50. Note, grabbing first instance since some have two sets, such
  #     as blood vs marrow samples. Also assuming that the EC50 value comes before the standard deviation
  #     column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  }
  col.auc <- which(str_detect(colnames(drug.data), coll("AUC", ignore_case = TRUE)))[1]

  # Gather the EC50 values if both they and the compound names exist in the file
  if (col.EC50 > 0 && col.compound > 0) {
    for (drug in drug.data[[col.compound]]){
      # Check to see if drug name is in compound list; if not, it's either:
      #     a) a synonym, or
      #     b) a removed compound
      #     Need to treat each category accordingly
      if (tolower(drug) %in% tolower(colnames(data.drug.ec50))) { # Drug name is in full compound dataset
        # Get index for drug compound in the full dataset
        data.drug.ec50.index <- which(tolower(colnames(data.drug.ec50)) == tolower(drug))
        # Get index for drug compound in the current patient's data
        # Note: only using value for first drug compound of the same name; some patient's have
        # duplicates, but one set will be blasts and one is something else
        drug.data.index <- which(drug.data[,col.compound] == drug)[1]
        # Get EC50 value if the column is called that; else get the IC50 value
        data.drug.ec50[patient.index, data.drug.ec50.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index, col.EC50])))
      } else { # Drug name is not in full compound dataset
        # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
        syn.index <- syns.row(drug)
        if (syn.index > 0) { # Must be a synonym drug
          # Preferred drug name in full set is the first name in synonym row
          # Get index for drug compound in the full dataset
          data.drug.ec50.index <- which(tolower(colnames(data.drug.ec50)) == tolower(synonyms[syn.index, 1]))
          drug.data.index <- which(drug.data[,col.compound] == drug)[1]
          # Get EC50 value if the column is called that; else get the IC50 value
          data.drug.ec50[patient.index, data.drug.ec50.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index,col.EC50])))
        } # else must be a removed drug and should do nothing with it
      }
    }
  } # else wasn't able to find the auc column and the patient data gathering was skipped

    # Gather the auc values if both they and the compound names exist in the file
  if (col.auc > 0 && col.compound > 0) {
    for (drug in drug.data[[col.compound]]){
      # Check to see if drug name is in compound list; if not, it's either:
      #     a) a synonym, or
      #     b) a removed compound
      #     Need to treat each category accordingly
      if (tolower(drug) %in% tolower(colnames(data.drug.auc))) { # Drug name is in full compound dataset
        # Get index for drug compound in the full dataset
        data.drug.auc.index <- which(tolower(colnames(data.drug.auc)) == tolower(drug))
        # Get index for drug compound in the current patient's data
        # Note: only using value for first drug compound of the same name; some patient's have
        # duplicates, but one set will be blasts and one is something else
        drug.data.index <- which(drug.data[,col.compound] == drug)[1]
        # Get EC50 value if the column is called that; else get the IC50 value
        data.drug.auc[patient.index, data.drug.auc.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index, col.auc])))
      } else { # Drug name is not in full compound dataset
        # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
        syn.index <- syns.row(drug)
        if (syn.index > 0) { # Must be a synonym drug
          # Preferred drug name in full set is the first name in synonym row
          # Get index for drug compound in the full dataset
          data.drug.auc.index <- which(tolower(colnames(data.drug.auc)) == tolower(synonyms[syn.index, 1]))
          drug.data.index <- which(drug.data[,col.compound] == drug)[1]
          # Get EC50 value if the column is called that; else get the IC50 value
          data.drug.auc[patient.index, data.drug.auc.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index,col.auc])))
        } # else must be a removed drug and should do nothing with it
      }
    }
  } # else wasn't able to find the auc column and the patient data gathering was skipped
}

# Go through all csv files and start collecting EC50 data
for (data.file in filenames.drug.csv) {

  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  patient.index <- which(rownames(data.drug.ec50) == patient)

  # Open file
  drug.data <- read.csv(fs::path(clean.drug.dir, data.file), check.names = FALSE, stringsAsFactors = FALSE)

  # Get column that has drug compounds, if exists; else give error
  if (length(which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))) > 0) {
    col.compound <- which(str_detect(colnames(drug.data), coll("compound", ignore_case = TRUE)))[1]
  }

  # Get column that has either EC50 or IC50. Note, grabbing first instance since some have two sets, such
  #     as blood vs marrow samples. Also assuming that the EC50 value comes before the standard deviation
  #     column that sometimes has EC50 in the name.
  if (length(which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("EC50", ignore_case = TRUE)))[1]
  } else if (length(which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))) > 0) {
    col.EC50 <- which(str_detect(colnames(drug.data), coll("IC50", ignore_case = TRUE)))[1]
  }

  # Gather the EC50 values if both they and the compound names exist in the file
  if (col.EC50 > 0 && col.compound > 0) {
    for (drug in drug.data[[col.compound]]){
      # Check to see if drug name is in compound list; if not, it's either:
      #     a) a synonym, or
      #     b) a removed compound
      #     Need to treat each category accordingly
      if (tolower(drug) %in% tolower(colnames(data.drug.ec50))) { # Drug name is in full compound dataset
        # Get index for drug compound in the full dataset
        data.drug.ec50.index <- which(tolower(colnames(data.drug.ec50)) == tolower(drug))
        # Get index for drug compound in the current patient's data
        # Note: only using value for first drug compound of the same name; some patient's have
        # duplicates, but one set will be blasts and one is something else
        drug.data.index <- which(drug.data[,col.compound] == drug)[1]
        # Get EC50 value if the column is called that; else get the IC50 value
        data.drug.ec50[patient.index, data.drug.ec50.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index, col.EC50])))
      } else { # Drug name is not in full compound dataset
        # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
        syn.index <- syns.row(drug)
        if (syn.index > 0) { # Must be a synonym drug
          # Preferred drug name in full set is the first name in synonym row
          # Get index for drug compound in the full dataset
          data.drug.ec50.index <- which(tolower(colnames(data.drug.ec50)) == tolower(synonyms[syn.index, 1]))
          drug.data.index <- which(drug.data[,col.compound] == drug)[1]
          # Get EC50 value if the column is called that; else get the IC50 value
          data.drug.ec50[patient.index, data.drug.ec50.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index,col.EC50])))
        } # else must be a removed drug and should do nothing with it
      }
    }
  } # else wasn't able to find the EC50 column and the patient data gathering was skipped
  
  col.auc <- which(str_detect(colnames(drug.data), coll("AUC", ignore_case = TRUE)))[1]
  # Gather the auc values if both they and the compound names exist in the file
  if (col.auc > 0 && col.compound > 0) {
    for (drug in drug.data[[col.compound]]){
      # Check to see if drug name is in compound list; if not, it's either:
      #     a) a synonym, or
      #     b) a removed compound
      #     Need to treat each category accordingly
      if (tolower(drug) %in% tolower(colnames(data.drug.auc))) { # Drug name is in full compound dataset
        # Get index for drug compound in the full dataset
        data.drug.auc.index <- which(tolower(colnames(data.drug.auc)) == tolower(drug))
        # Get index for drug compound in the current patient's data
        # Note: only using value for first drug compound of the same name; some patient's have
        # duplicates, but one set will be blasts and one is something else
        drug.data.index <- which(drug.data[,col.compound] == drug)[1]
        # Get EC50 value if the column is called that; else get the IC50 value
        data.drug.auc[patient.index, data.drug.auc.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index, col.auc])))
      } else { # Drug name is not in full compound dataset
        # Try to get the row index for the drug from synonym file; will be 0 if no synonyms
        syn.index <- syns.row(drug)
        if (syn.index > 0) { # Must be a synonym drug
          # Preferred drug name in full set is the first name in synonym row
          # Get index for drug compound in the full dataset
          data.drug.auc.index <- which(tolower(colnames(data.drug.auc)) == tolower(synonyms[syn.index, 1]))
          drug.data.index <- which(drug.data[,col.compound] == drug)[1]
          # Get EC50 value if the column is called that; else get the IC50 value
          data.drug.auc[patient.index, data.drug.auc.index] <- suppressWarnings(suppressMessages(as.numeric(drug.data[drug.data.index,col.auc])))
        } # else must be a removed drug and should do nothing with it
      }
    }
  } # else wasn't able to find the auc column and the patient data gathering was skipped
}

##### Drug Sensitivity Data Thresholds & Transform #####

# Set thresholds on EC50 values
# Any value > active.maximum is set to active.maximum/1e-3 to make large enough to be different from barely active drugs
# Any value < active.minimum is set to active.minimum
# Upper bound threshold to set inactive values to
drug.upper.threshold <- opts$active.maximum/1e-3

drugs <- data.drug.ec50

# Set all threshold values, if needed
for (patient in rownames(data.drug.ec50)) {
  for (drug in colnames(data.drug.ec50)) {
    # Make sure not NA
    if (!is.na(data.drug.ec50[patient, drug])) {
      # Check if value is higher than active threshold and set to upper threshold value, if so
      if (data.drug.ec50[patient, drug] > opts$active.maximum) {
        data.drug.ec50[patient, drug] <- drug.upper.threshold
      } else if (data.drug.ec50[patient, drug] < opts$active.minimum) {
        data.drug.ec50[patient, drug] <- opts$active.minimum
      }
    } # else leave as NA
  }
}

# -log(EC50) of all EC50 values
data.drug.ec50 <- apply(data.drug.ec50, 2, function(x) {-log(x)})

##### Variant Data Gathering #####

# Get list of the files in the directory
filenames.var <- list.files(clean.var.dir, pattern = "*.xlsx")

# Patient list from variant cleaning
patients.var <- as.numeric(readLines(fs::path(opts$var.dir, "Patients_Variants.txt")))
# Gene list from Genes_Remove.csv
genes.remove <- read.csv(fs::path(opts$var.dir, "Genes_Remove.csv"), header=TRUE, na.strings = c("", " ", NA), stringsAsFactors = FALSE, check.names = FALSE)
genes <- genes.remove$Gene[which(genes.remove$Include == 1)]

# Dataframe for the collected variant data
data.var <- data.frame(matrix(nrow=length(patients.var), ncol=length(genes)))
rownames(data.var) <- patients.var
colnames(data.var) <- genes
# Initialize all values to 0
data.var[is.na(data.var)] <- 0

# Go through all files again and start collecting the variant data
for (data.file in filenames.var) {

  # Get names of sheets in Excel workbook
  sheets <- excel_sheets(fs::path(clean.var.dir,data.file))

  sheet.miss <- which(str_detect(sheets, coll("missense", ignore_case = TRUE)))
  sheet.indel <- which(str_detect(sheets, coll("indel", ignore_case = TRUE)))

  # Get patient ID number from file name
  patient <- as.numeric(str_extract(data.file, "[0-9]+"))
  patient.index <- which(rownames(data.var) == patient)

  # Open file to missense info
  var.data <- read_excel(fs::path(clean.var.dir, data.file), sheet=sheets[sheet.miss])

  # Filter based on startup arguments
  # SIFT filtering
  sift.ind <- vector()
  if (opts$sift.deleterious) {
    # Get indices for deleterious SIFT
    sift.ind <- c(sift.ind, grep("delet", var.data$SIFT))
  } # Else not chosen and nothing added to sift.ind
  if (opts$sift.tolerated) {
    # Get indices for tolerated SIFT
    sift.ind <- c(sift.ind, grep("toler", var.data$SIFT))
  } # Else not chosen and nothing added to sift.ind
  # Only keep indices gathered in filters above
  # To be safe, check that sift.ind is not empty; this could happen if neither option is checked
  if (length(sift.ind) > 0) {
    var.data <- var.data[sift.ind,]
  } # else keep all data

  # Polyphen filtering
  poly.ind <- vector()
  if (opts$poly.prob.damaging) {
    # Get indices for probably damaging polyphen
    poly.ind <- c(poly.ind, grep("prob", var.data$Polyphen))
  } # Else not chosen and nothing added to poly.ind
  if (opts$poly.poss.damaging) {
    # Get indices for possibly damaging polyphen
    poly.ind <- c(poly.ind, grep("poss", var.data$Polyphen))
  } # Else not chosen and nothing added to poly.ind
  if (opts$poly.benign) {
    # Get indices for benign polyphen
    poly.ind <- c(poly.ind, grep("benign", var.data$Polyphen))
  }
  # Only keep indices gathered in filters above
  # To be safe, check that poly.ind is not empty; this could happen if no options are checked
  if (length(poly.ind) > 0) {
    var.data <- var.data[poly.ind,]
  } # else keep all data

  # Remove missense with freq < opts$Frequency (default is 2.5%)
  # Note that currently expect value in %, not decimal so / by 100
  var.data <- subset(var.data, var.data$Frequency > (opts$frequency/100) )

  # Count number of variations in each gene and store
  for (gene in var.data$Gene){
    # Get index for gene in the full dataset
    data.var.index <- which(colnames(data.var) == gene)
    # Add 1 to current value
    data.var[patient.index, data.var.index] <- data.var[patient.index, data.var.index] + 1
  }

  # Open file to indel info
  var.data <- read_excel(fs::path(clean.var.dir, data.file), sheet=sheets[sheet.indel])
  # Remove indel with freq < opts$Frequency (default is 2.5%)
  # Note that currently expect value in %, not decimal so / by 100
  var.data <- subset(var.data, var.data$Frequency > (opts$frequency/100) )

  # Count number of variations in each gene and store
  for (gene in var.data$Gene){
    # Get index for gene in the full dataset
    data.var.index <- which(colnames(data.var) == gene)
    # Add 1 to current value
    data.var[patient.index, data.var.index] <- data.var[patient.index, data.var.index] + 1
  }
}

##### Subsetting data by which patients have both drug and variant data
patient.subset <- intersect(rownames(data.drug.ec50), rownames(data.var))
# Assuming that the subset has at least one, but put a condition just in case.
if (length(patient.subset) > 0) {
  data.var <- data.var[patient.subset,]
  data.drug.ec50 <- data.drug.ec50[patient.subset,]
} else {
  stop("There appears to be no shared patients in the dataset.")
}

##### Subsetting data by patient list #####

# Need to make subsetting file if it doesn't exist
#   Subsetting file is for getting the patient subset the user wants
# Check to see if Patient_Subset.csv has been created, yet
#    If no, get patient list that has both data sets, create file, and ask user to edit subset
if (!fs::file_exists(fs::path(opts$output.dir, "Patient_Subset.csv"))) {
  # File doesn't exist yet; create it
  # Create a dataframe for patient list and whether to include
  patients.include <- data.frame(matrix(nrow=length(patient.subset), ncol=2))
  colnames(patients.include) <- c("Patient", "Include")
  patients.include$Patient <- patient.subset
  # Default is to include all patients
  patients.include$Include <- 1
  # Write to file
  write.csv(patients.include, fs::path(opts$output.dir, "Patient_Subset.csv"), row.names = FALSE)
} else {
  # Create a dataframe for patient list and whether to include
  patients.include <- data.frame(matrix(nrow=length(patient.subset), ncol=2), check.names = FALSE)
  colnames(patients.include) <- c("Patient", "Include")
  patients.include$Patient <- patient.subset
  # Default is to include all patients
  patients.include$Include <- 1
  # Open the file
  patient.subset.file <- read.csv(fs::path(opts$output.dir, "Patient_Subset.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  patient.subset.exclude <- patient.subset.file[which(patient.subset.file$Include == 0), ]
  for (patient in patient.subset.exclude$Patient) {
    index <- which(rownames(data.var) == patient)
    if (length(index) > 0) {
      patients.include$Include[index] <- 0
    }
  }
  # Write to file
  write.csv(patients.include, fs::path(opts$output.dir, "Patient_Subset.csv"), row.names = FALSE)
}

##### Subset data by drug compound list #####

# Only need the subset of compounds that the user wants
# Check to see if Compound_Subset.csv has been created, yet
#    If no, get compound list that is not NA for the subset of patients chosen
if (!fs::file_exists(fs::path(opts$output.dir, "Compounds_Remove.csv"))) {
  # File doesn't exist, yet
  # Create a dataframe for compound list and whether to include
  compounds.include <- data.frame(matrix(nrow=dim(data.drug.ec50)[2], ncol=2))
  colnames(compounds.include) <- c("Compound", "Include")
  compounds.include$`Compound` <- colnames(data.drug.ec50)
  # Default is to include all compounds
  compounds.include$Include <- 1
  # Write to file
  write.csv(compounds.include, fs::path(opts$output.dir, "Compound_Subset.csv"), row.names = FALSE)
} else {
  # Create a dataframe for compound list and whether to include
  compounds.include <- data.frame(matrix(nrow=dim(data.drug.ec50)[2], ncol=2))
  colnames(compounds.include) <- c("Compound", "Include")
  compounds.include$`Compound` <- colnames(data.drug.ec50)
  # Default is to include all compounds
  compounds.include$Include <- 1
  compound.subset.file <- read.csv(fs::path(opts$output.dir, "Compound_Subset.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  compound.subset.exclude <- compound.subset.file[which(compound.subset.file$Include == 0), ]
  for (compound in compound.subset.exclude$`Compound`) {
    index <- which(colnames(data.var) == compound)
    if (length(index) > 0) {
      compounds.include$Include[index] <- 0
    }
  }
  # Write to file
  write.csv(compounds.include, fs::path(opts$output.dir, "Compound_Subset.csv"), row.names = FALSE)
}

##### Write to files #####

# Output filtered, subset, etc tables to files
write.csv(data.drug.ec50, file=fs::path(opts$output.dir, "PatientsVsDrugs.csv"))
write.csv(data.drug.auc, file=fs::path(opts$output.dir, "PatientsVsDrugs_auc.csv"))
write.csv(data.var, file=fs::path(opts$output.dir, "PatientsVsVariants.csv"))

# Write config file
write.csv(list(opts), fs::path(opts$output.dir, "config.csv"))
