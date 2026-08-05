### Made by Zach Burcham for correcting batch effect caused by different studies in the BSFL gut microbiome meta-analysis

###
#
# Run this program after running ConQuR-Default.R to generate the PERMANOVA stats and file
#
###

############################################################################################
#
# PERMANOVA R2
#
############################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

# PERMANOVA Calculation has problem running after the previous steps unless you clear and reload your env
invisible(gc()) # clear memory
rm(list = ls()) #clear environment

# package names
packages <- c("devtools","coda.base", "doParallel", "tictoc", "openxlsx", 
              "ggplot2", "dplyr", "reshape2", "readr", 
              "qiime2R", "phyloseq", "ape","randomForest",
              "caret","MLmetrics","multcompView","optparse")

# Function to check if each package is installed and install if necessary
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  # Load the package
  library(pkg, character.only = TRUE)
}
rm(packages,pkg)

# Install ConQuR with the ADONIS commands changed to ADONIS2

install.packages("/lustre/isaac24/scratch/rsaho/tools/ConQuR", repos = NULL, type = "source")
library(ConQuR)



# Define command-line arguments using optparse
option_list <- list(
  make_option(c("--output_dir"), type = "character", default = NULL, 
              help = "Path to the output directory where results will be stored", metavar = "character")
)

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Check if necessary files are provided
if (is.null(opt$output_dir)) {
  stop("You must provide input: --output_dir")
}

# Create the base output directory if it doesn't exist
output_dir <- opt$output_dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message(paste("Output directory", output_dir, "created."))
} else {
  message(paste("Output directory", output_dir, "already exists."))
}

# Load previous data
load(file = file.path(file.path("batch_correction", "checkpoints", "all_studies_corrections.RData")))

# Run PERMANOVA analysis and store results
cat("Beginning PERMANOVA R2 Calculation\n")
tic("Completed PERMANOVA R2 Calculation")

# Create empty lists to store permanova r2 results
permanova_standard_tabcount_results <- list()
permanova_sqrt_dist_tabcount_results <- list()
permanova_add_tabcount_results <- list()

permanova_standard_tabrel_results <- list()
permanova_sqrt_dist_tabrel_results <- list()
permanova_add_tabrel_results <- list()

# Get a list of corrections
corr_list <- ls(pattern = "^taxa_corrected_ref_.*")
pen_corr_list <- ls(pattern = "^taxa_pen_corrected_ref_.*")

# Calculate permanoca R2 for uncorrected
cat("Calculating PERMANOVA R2 for Original Data\n")
permanova_before = suppressWarnings(PERMANOVA_R2(TAX=taxa, batchid=batchid, covariates=covar, key_index=1))
cat("Calculation Complete\n")

# Loop through corrections and calculate permanova R2
for (study in batchid_list) {
  cat(paste("Calculating PERMANOVA R2 for Corrected Data using Ref. Study", study, "\n"))
  corr_id <- paste0("taxa_corrected_ref_",study)
  pen_corr_id <- paste0("taxa_pen_corrected_ref_",study)
  corr_matrix <- get(corr_id)
  pen_corr_matrix <- get(pen_corr_id)
  permanova_after_default = suppressWarnings(PERMANOVA_R2(TAX=corr_matrix, batchid=batchid, covariates=covar, key_index=1))
  permanova_after_penalized = suppressWarnings(PERMANOVA_R2(TAX=pen_corr_matrix, batchid=batchid, covariates=covar, key_index=1))
  # Extract and combine tabcount PERMANOVA results
  permanova_standard_tabcount_results[[study]] <- data.frame(
    Study = study,
    Correction = c("None", "Default", "Penalized"),
    Batch_Standard = c(permanova_before$tab_count["batch", "standard"],
                       permanova_after_default$tab_count["batch", "standard"],
                       permanova_after_penalized$tab_count["batch", "standard"]),
    Key_Standard = c(permanova_before$tab_count["key", "standard"],
                     permanova_after_default$tab_count["key", "standard"],
                     permanova_after_penalized$tab_count["key", "standard"])
  )

  permanova_sqrt_dist_tabcount_results[[study]] <- data.frame(
    Study = study,
    Correction = c("None", "Default", "Penalized"),
    Batch_Sqrt_Dist = c(permanova_before$tab_count["batch", "sqrt.dist=T"],
                        permanova_after_default$tab_count["batch", "sqrt.dist=T"],
                        permanova_after_penalized$tab_count["batch", "sqrt.dist=T"]),
    Key_Sqrt_Dist = c(permanova_before$tab_count["key", "sqrt.dist=T"],
                      permanova_after_default$tab_count["key", "sqrt.dist=T"],
                      permanova_after_penalized$tab_count["key", "sqrt.dist=T"])
  )

  permanova_add_tabcount_results[[study]] <- data.frame(
    Study = study,
    Correction = c("None", "Default", "Penalized"),
    Batch_Add = c(permanova_before$tab_count["batch", "add=T"],
                  permanova_after_default$tab_count["batch", "add=T"],
                  permanova_after_penalized$tab_count["batch", "add=T"]),
    Key_Add = c(permanova_before$tab_count["key", "add=T"],
                permanova_after_default$tab_count["key", "add=T"],
                permanova_after_penalized$tab_count["key", "add=T"])
  )

  # Extract and combine tabrel PERMANOVA results
  permanova_standard_tabrel_results[[study]] <- data.frame(
    Study = study,
    Correction = c("None", "Default", "Penalized"),
    Batch_Standard = c(permanova_before$tab_rel["batch", "standard"],
                       permanova_after_default$tab_rel["batch", "standard"],
                       permanova_after_penalized$tab_rel["batch", "standard"]),
    Key_Standard = c(permanova_before$tab_rel["key", "standard"],
                     permanova_after_default$tab_rel["key", "standard"],
                     permanova_after_penalized$tab_rel["key", "standard"])
  )

  permanova_sqrt_dist_tabrel_results[[study]] <- data.frame(
    Study = study,
    Correction = c("None", "Default", "Penalized"),
    Batch_Sqrt_Dist = c(permanova_before$tab_rel["batch", "sqrt.dist=T"],
                        permanova_after_default$tab_rel["batch", "sqrt.dist=T"],
                        permanova_after_penalized$tab_rel["batch", "sqrt.dist=T"]),
    Key_Sqrt_Dist = c(permanova_before$tab_rel["key", "sqrt.dist=T"],
                      permanova_after_default$tab_rel["key", "sqrt.dist=T"],
                      permanova_after_penalized$tab_rel["key", "sqrt.dist=T"])
  )

  permanova_add_tabrel_results[[study]] <- data.frame(
    Study = study,
    Correction = c("None", "Default", "Penalized"),
    Batch_Add = c(permanova_before$tab_rel["batch", "add=T"],
                  permanova_after_default$tab_rel["batch", "add=T"],
                  permanova_after_penalized$tab_rel["batch", "add=T"]),
    Key_Add = c(permanova_before$tab_rel["key", "add=T"],
                permanova_after_default$tab_rel["key", "add=T"],
                permanova_after_penalized$tab_rel["key", "add=T"])
  )
  cat("Calculation Complete\n")
}
toc()
save(list = ls(), file = file.path(output_dir, "checkpoints", "default_permanova_results.RData"))


cat("Compiling PERMANOVA R2 Results\n")
tic("Compiled PERMANOVA R2 Results: Saved to permanova/final_permanova_results.xlsx")

# Combine all tabcount results into single data frames
permanova_standard_tabcount_combined <- do.call(rbind, permanova_standard_tabcount_results)
permanova_sqrt_dist_tabcount_combined <- do.call(rbind, permanova_sqrt_dist_tabcount_results)
permanova_add_tabcount_combined <- do.call(rbind, permanova_add_tabcount_results)

permanova_standard_tabrel_combined <- do.call(rbind, permanova_standard_tabrel_results)
permanova_sqrt_dist_tabrel_combined <- do.call(rbind, permanova_sqrt_dist_tabrel_results)
permanova_add_tabrel_combined <- do.call(rbind, permanova_add_tabrel_results)

# Merge into final data frames and save
final_tabcount_permanova_combined <- merge(merge(permanova_standard_tabcount_combined, permanova_sqrt_dist_tabcount_combined, by=c("Study", "Correction")), permanova_add_tabcount_combined, by=c("Study", "Correction"))
final_tabrel_permanova_combined <- merge(merge(permanova_standard_tabrel_combined, permanova_sqrt_dist_tabrel_combined, by=c("Study", "Correction")), permanova_add_tabrel_combined, by=c("Study", "Correction"))

# Calculate the batch/key ratio for tabcount results
final_tabcount_permanova_combined$Batch_Key_Ratio_Standard <- (final_tabcount_permanova_combined$`Batch_Standard` / final_tabcount_permanova_combined$`Key_Standard`)
final_tabcount_permanova_combined$Batch_Key_Ratio_Sqrt_Dist <- (final_tabcount_permanova_combined$`Batch_Sqrt_Dist` / final_tabcount_permanova_combined$`Key_Sqrt_Dist`)
final_tabcount_permanova_combined$Batch_Key_Ratio_Add <-  (final_tabcount_permanova_combined$`Batch_Add` / final_tabcount_permanova_combined$`Key_Add`)

# Calculate the batch/key ratio  for tabrel results
final_tabrel_permanova_combined$Batch_Key_Ratio_Standard <- (final_tabrel_permanova_combined$`Batch_Standard` / final_tabrel_permanova_combined$`Key_Standard`)
final_tabrel_permanova_combined$Batch_Key_Ratio_Sqrt_Dist <- (final_tabrel_permanova_combined$`Batch_Sqrt_Dist` / final_tabrel_permanova_combined$`Key_Sqrt_Dist`)
final_tabrel_permanova_combined$Batch_Key_Ratio_Add <- (final_tabrel_permanova_combined$`Batch_Add` / final_tabrel_permanova_combined$`Key_Add`)

# Save permanova results in a new excel workbook
wb <- createWorkbook()

# Add worksheets for each data frame
addWorksheet(wb, "TabCount_PERMANOVA_R2")
addWorksheet(wb, "TabRel_PERMANOVA_R2")

# Write the data frames to their respective worksheets
writeData(wb, sheet = "TabCount_PERMANOVA_R2", final_tabcount_permanova_combined)
writeData(wb, sheet = "TabRel_PERMANOVA_R2", final_tabrel_permanova_combined)

# Save the workbook as an Excel file
saveWorkbook(wb, file.path(output_dir, "permanova", "final_permanova_results.xlsx"), overwrite = TRUE)
toc()

save(list = ls(), file = file.path(output_dir, "checkpoints", "default_permanova_program_complete.RData"))

cat("***\nProgram Complete\n***")
