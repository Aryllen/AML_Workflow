# Workflow for Discovery of Gene Variants that Drive Drug Sensitivity in AML Patients

This workflow is for [Biodepot Workflow Builder (BwB)](https://github.com/BioDepot/BioDepot-workflow-builder). The workflow requires genetic variant data from [MyAML](https://catalog.invivoscribe.com/product/myaml-194-targeted-ngs-gene-panel-clia/) and drug sensitivity lab results. The output is a Jupyter Notebook analysis, comparing drug sensitivity of patients with a given genetic variant against patients without the genetic variant. Given that BwB works completely within Docker containers, not relying on the host system for programming language or package requirements, this workflow produces repeatable results for the same dataset on any machine that can run the workflow within BwB.

## Setup Requirements

### Docker

Docker is a requirement of BwB. Information on how to install the community version of Docker can be found at the [Docker website](https://docs.docker.com/v17.12/install/#desktop) or within the [BwB Readme](https://github.com/BioDepot/BioDepot-workflow-builder).

### Biodepot Workflow Builder (BwB)

Visit the [BwB Readme](https://github.com/BioDepot/BioDepot-workflow-builder) for how to run BwB in Docker.

### AML_Workflow

Of course, the AML_Workflow is also required. [Clone](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github) this repo any local directory that is able to be seen from within the BwB docker container.

## Data Requirements

There are some expectations regarding the input data for the workflow. Please ensure the data matches the guidelines below. With the exception of file types, most of these requirements are checked for in the cleaning section of the workflow. However, verifying that the data is properly formatted ahead of time will make for less issues for the workflow to find.

**NOTE:** All docker images in the workflow are built on a flavor of Linux. As such, it is expected that there are no spaces in any directory or file name that is being accessed. This includes any directory along the file path.

### Variant Data

At this time, the variant data is expected to be from MyAML and in the form of an Excel spreadsheet. There should be one file per patient and a numeric patient ID should be in the filename. Two sheets, in particular, are expected: one with the word "missense" in the name, and one with the word "indel" in the name. The indel worksheet requires the following columns: Gene, Frequency. The missense worksheet requires the following columns: Gene, Frequency, Polyphen, SIFT, Sample. The Sample file should have a string with the numeric patient ID near the beginning and separated by some type of punctuation from any other numeric value.

### Drug Sensitivity Data

The drug sensitivity data can be in either an Excel spreadsheet or csv file. There should be one file for each patient and a numeric patient ID should be in the filename. For Excel spreadsheets, the data should be on a sheet with one of the following phrases in the name: 
- fitted parameters static
- fitparameters
- combined
- blast param statsort
- parameters static
- 4database

Drug sensitivity data should have the following required columns: Compound, either Sample or File, and either EC50 or IC50. The Sample or File column should have the numeric patient ID alone or a string with the patient ID being the first set of numbers that appear in the string (any other numeric value, if present, in the string should be separated by some form of punctuation from the ID).

## Running the Workflow

Opening the workflow in BwB:

1. Start the BwB container. Instructions for running the container can be found in the [BwB Readme](https://github.com/BioDepot/BioDepot-workflow-builder).
2. Select 'Load Workflow' from the 'File' dropdown in the main menu.
3. Navigate to the cloned AML_Workflow directory and click 'Choose'.

The workflow should now be visible in the BwB Workspace. Three widgets require input before they start, Variant Initial Cleaning, Drug Sensitivity Initial Cleaning, and Initial Analysis. Double-click on each of these to input the needed directories and, in the case of Initial Analysis, any optional parameters desired. More detailed information on these widgets can be found in the Widgets section of this Readme. Once the required fields are completed, double-click on Variant Initial Cleaning and Drug Sensitivity Initial Cleaning and select 'Start'. This will begin running the workflow. Read more about what should be expected during the workflow for each branch in the Widgets section below.

## Widgets

### Variant Branch

#### Initial Variant Cleaning

Input: directory with individual patient variant files (variant_dir)
Output: 
- Always:
  - variant_dir/output (an output directory within given variant directory)
- If data has 1+ errors:
  - variant_dir/output/Variant_File_Errors.txt
- If data has no errors:
  - variant_dir/output/Genes_Remove.csv*
  - variant_dir/output/\_Clean
    - This directory should never be touched. Clean variant files will be copied to this directory and used in the Initial Analysis widget. Putting unclean files in this directory may result in the workflow failing.
  * if file does not exist, yet

The Initial Variant Cleaning (IVC) widget checks each file within given variant directory to ensure the following:
- Sample column exists,
- numeric patient ID in filename matches numeric patient ID from Sample column,
- a sheet with the word "missense" in the name exists,
- a sheet with the word "indel" in the name exists.

If any of the files have an error, the IVC widget will output a text file (Variant_File_Errors.txt) with the errors found. The text file can be found within the new Output directory within the original variant directory. Additionally, this will signal the error branch of the workflow, which pops up the text file in gedit via the Variant File Errors widget and stops the workflow from continuing further. The errors found must be rectified before continuing. Once the errors are fixed, restart the IVC widget to verify cleanliness.

#### Show Errors If Variant Data Unclean/Verify Variant Data Clean

These two widgets act as gatekeepers for the workflow. If the variant data is found to have errors in the IVC widget, then Verify Variant Data Clean (VDC) will not allow the workflow to continue towards the analysis branch. Instead, Show Errors If Variant Data Unclean (VDU) will trigger the next widget, Variant File Errors, to show the list of errors found in the files. Conversely, if the variant data has no errors, the VDU widget will halt the error branch and the VDC widget will allow the workflow to continue onto the Edit Genes for Removal widget.

#### Variant File Errors

This widget opens gedit for viewing the Variant_File_Errors.txt file. This file can also be located within the given variant directory in the created output directory.

#### Edit Genes for Removal

This widget opens the Genes_Remove.csv file that either existed from a previous run of the workflow or was created in the Initial Variant Cleaning widget. There are two columns, Gene and Include. The Include column is binary, in that 1 means to include the gene in further analysis and 0 means to not include the gene. This is the time to remove any genes that simply should not be included in the analysis by setting Include for that row to 0. After saving (see important note below), exit gnumeric to set the variant trigger for the Initial Analysis widget.

**Important:** The file must be saved as a csv, which is a little complicated in gnumeric. To do this, click on 'Data' in the main menu toolbar, then 'Export Data' and 'Export as CSV'. Click 'Save' and select 'Okay' to overwrite the file. Do not change the filename. Exit gnumeric and select 'discard changes' when prompted.

### Drug Sensitivity Branch

#### Initial Drug Sensitivity Cleaning

Input: directory with individual patient drug sensitivity files (drug_dir)
Output: 
- Always:
  - drug_dir/output (an output directory within drug sensitivity directory)
- If data has 1+ errors:
  - drug_dir/output/Drug_Sensitivity_File_Errors.txt
- If data has no errors:
  - drug_dir/output/Compounds_Remove.csv*
  - drug_dir/output/Compound_Synonyms.csv*
  - drug_dir/output/\_Clean
    - This directory should never be touched. Clean drug sensitivity files will be copied to this directory and used in the Initial Analysis widget. Putting unclean files in this directory may result in the workflow failing.
  * if files do not exist, yet


The Initial Drug Sensitivity Cleaning (IDSC) widget checks each file within the given drug sensitivity directory to ensure the following:
- a Sample or File column exists,
- a column called either EC50 or IC50 exists,
- a Compound column exists,
- numeric patient ID in filename matches the numeric patient ID found in either the Sample or File column,
- a sheet with one of the accepted phrases (see Drug Sensitivity Data Requirements above) in the name exists,
- data in the EC50 or IC50 column is numeric.

If any of the files have an error, the IDSC widget will output a text file (Drug_Sensitivity_File_Errors.txt) with the errors found. The text file can be found within the new Output directory within the original drug sensitivity directory. Additionally, this will signal the error branch of the workflow, which pops up the text file in gedit via the Variant File Errors widget and stops the workflow from continuing further. The errors found must be rectified before continuing. Once the errors are fixed, restart the IVC widget to verify cleanliness.

#### Drug Sensitivity File Errors

This widget opens gedit for viewing the Drug_Sensitivity_File_Errors.txt file. This file can also be located on your local drive that was mapped to BwB under the created output directory within the given drug sensitivity directory.

#### Edit Drug Synonym List

The Edit Drug Synonym List exists to create a dictionary of preferred drug names and synonym names that may/do show up in the drug sensitivity data. This is to help with data entry errors and changes in drug names (e.g. brand name versus generic name). If the synonym names are not accounted for in this file, then each drug name is considered a separate, unique drug in the analysis.

The Edit Drug Synonym List widget opens gnumeric with the Compound_Synonyms.csv file. The columns are not named in this file. The first column is the preferred name for a drug, while the synonym names go in the next columns for the row, one synonym per column. The preferred name is case-sensitive, but the synonym names are not. 

For example:
| True name | Synonym 1 | Synonym 2 |
| ------------- | -------------------------- | --------- |
| larotrectinib |	Larotrectinib (LOXO-101) |	loxo-101 |
| Venetoclax    |	Venetoclax (ABT-199)     |	abt-199  |

In the table above, larotrectinib and Venetoclax would be the preferred drug names and would appear this way in the later analysis. When gathering the EC50/IC50 values from the drug sensitivity files, any files with Larotrectinib (LOXO-101) or loxo-101 would be stored under the preferred name larotrectinib. Similarly for Venetoclax (ABT-199) and abt-199; the values for these compounds would be stored under Venetoclax.

The easiest way to see if there are drugs that should have preferred names and synonym names accounted for in the dataset is to navigate to the Compounds_Remove.csv created in the Initial Drug Sensitivity Cleaning widget. This file can be found within the created output directory of the given drug sensitivity directory. Open this file and scan the compounds in the Compound column for names that refer to the same drug. The file is in alphabetical order making data entry mistakes easy to find. However, finding generic names that should be under the brand name (or vice versa) requires knowledge of both drug names.

Once the Compound_Synonyms.csv file has been updated (see important note below), exit gnumeric to trigger the next widget in the workflow, Clean Drugs of Synonyms.

**Important:** The file must be saved as a csv, which is a little complicated in gnumeric. To do this, click on 'Data' in the main menu toolbar, then 'Export Data' and 'Export as CSV'. Click 'Save' and select 'Okay' to overwrite the file. Do not change the filename. Exit gnumeric and select 'discard changes' when prompted.

#### Clean Drugs of Synonyms

The Clean Drugs of Synonyms widget cleans out the synonyms from the Compounds_Remove.csv based on the Compound_Synonyms.csv dictionary created in Edit Compound Synonyms List widget.

#### Edit Drugs for Removal

This widget opens the Compounds_Remove.csv file that was created in the Initial Drug Sensitivity Cleaning widget. There are two columns, Compound and Include. The Include column is binary, in that 1 means to include the drug in further analysis and 0 means to not include the drug. This is the time to remove any drugs that simply should not be included in the analysis by setting Include for that row to 0. After saving (see important note below), exit gnumeric to set the drug sensitivity trigger for the Initial Analysis widget.

**Important:** The file must be saved as a csv, which is a little complicated in gnumeric. To do this, click on 'Data' in the main menu toolbar, then 'Export Data' and 'Export as CSV'. Click 'Save' and select 'Okay' to overwrite the file. Do not change the filename. Exit gnumeric and select 'discard changes' when prompted.

### Analysis Branch

#### Initial Analysis

Input: desired Output directory for analysis data (analysis_dir)
Output:
  - analysis_dir/PatientsVsVariants.csv
  - analysis_dir/PatientsVsDrugs.csv
  - analysis_dir/Patient_Subset.csv
  - analysis_dir/Compound_Subset.csv
  - analysis_dir/config.yml


The Initial Analysis (IA) widget gathers the variant and drug sensitivity data into two tables and outputs them to csv files: PatientsVsVariants.csv and PatientsVsDrugs.csv. The only required input to the widget is the output directory, analysis_dir, which is where files created by this widget will be written to, including the two csv files mentioned. There are several optional inputs that not only affect which data is collected, but also what happens to the data.

Due to the need for the Analysis Notebook to have patients that have both drug sensitivity and variant data, an initial cleaning removes all patients that do not have one or the other dataset.

#### Edit Patient Subset

The Edit Patient Subset widget opens the Patient_Subset.csv created by the previous widget, Initial Analysis. The csv has two columns, the Patients that have both drug and variant data and Include, a binary value for inclusion. Any patients who's Include value is changed to 0 will be removed from the PatientsVsVariants.csv by the Re-subset Patients & Compounds widget. This allows for getting analysis/testing on a limited subset of patients. Rerunning the previous widget, Initial Analysis, will get the full dataset, again, and a new subset can be chosen via the Edit Patient Subset widget. Once the patients are selected or removed, save the csv, and exit the widget to trigger the next widget.

**Important:** The file must be saved as a csv, which is a little complicated in gnumeric. To do this, click on 'Data' in the main menu toolbar, then 'Export Data' and 'Export as CSV'. Click 'Save' and select 'Okay' to overwrite the file. Do not change the filename. Exit gnumeric and select 'discard changes' when prompted.

#### Edit Drug Subset

The Edit Drug Subset widget opens the Compound_Subset.csv created by the previous widget, Initial Analysis. The csv has two columns, Compounds that at least one patient in the dataset had data for and Include, a binary value for inclusion. Any compounds with an Include value of 0 will be removed from the PatientsVsDrugs.csv by the Re-subset Patients & Compounds widget. This allows for getting analysis/testing on a limited subset of drugs. Rerunning the previous widget, Initial Analysis, will get the full dataset, again, and a new subset can be chosen via the Edit Drug Subset widget. Once the drugs are selected or removed, save the csv, and exit the widget to trigger the next widget.

**Important:** The file must be saved as a csv, which is a little complicated in gnumeric. To do this, click on 'Data' in the main menu toolbar, then 'Export Data' and 'Export as CSV'. Click 'Save' and select 'Okay' to overwrite the file. Do not change the filename. Exit gnumeric and select 'discard changes' when prompted.

#### Re-subset Patients & Compounds

Re-subset Patients & Compounds removes the compounds or patients with an Include value of 0. The widget will also remove compounds that have all NA values; this would happen in the case that patients removed from selection were the only ones with data for that drug.

#### Analysis Notebook

This widget opens a fully function Jupyter Notebook with the analysis.

#### Note

Multivariate and univariate model analysis is currently being done in a script. This script can be found in the Dockers/Analysis/multi_uni_analysis.R file.