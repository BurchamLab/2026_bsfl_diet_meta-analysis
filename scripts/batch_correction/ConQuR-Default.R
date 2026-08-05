### Made by Zach Burcham for correcting batch effect caused by different studies in the BSFL gut microbiome meta-analysis

invisible(gc()) # clear memory
rm(list = ls()) #clear environment

# check if these are installed and install if not, then load
# has to be run first and seperately because they need to be loaded to install the github packages
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# if (!requireNamespace("coda.base", quietly = TRUE)) install.packages("coda.base")
# if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
library(devtools)
library(coda.base)
library(RColorBrewer)

## install these seperately since they go through github, first ConQur is the original and the 2nd is a fork correcting the parallelization
# if (!requireNamespace("ConQuR", quietly = TRUE)) {
#   devtools::install_github("ivartb/ConQuR_par", build_vignettes = TRUE, force=TRUE)
# }

# if (!requireNamespace("qiime2R", quietly = TRUE)) {
#   install.packages("BiocManager")
#   BiocManager::install(c("Biostrings", "rhdf5", "phyloseq", "S4Vectors", "TreeSummarizedExperiment"))
#   devtools::install_github("jbisanz/qiime2R")
# }

# package names
packages <- c("ConQuR", "doParallel", "tictoc", "openxlsx", 
              "ggplot2", "dplyr", "reshape2", "readr", 
              "qiime2R", "phyloseq", "ape","randomForest",
              "caret","MLmetrics","multcompView","optparse", "RColorBrewer")

# Function to check if each package is installed and install if necessary
for (pkg in packages) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)
#   }
  # Load the package
  library(pkg, character.only = TRUE)
}
rm(packages,pkg)

# Define command-line arguments using optparse
option_list <- list(
  make_option(c("--tree"), type = "character", default = NULL, 
              help = "Path to the phylogenetic tree file (.qza)", metavar = "character"),
  make_option(c("--table"), type = "character", default = NULL, 
              help = "Path to the taxa table file (.qza)", metavar = "character"),
  make_option(c("--metadata"), type = "character", default = NULL, 
              help = "Path to the R-readable metadata file (tab delim)", metavar = "character"),
  make_option(c("--output_dir"), type = "character", default = NULL, 
              help = "Path to the output directory where results will be stored", metavar = "character")
)

# Parse the command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Check if necessary files are provided
if (is.null(opt$tree) || is.null(opt$table) || is.null(opt$metadata) || is.null(opt$output_dir)) {
  stop("You must provide all inputs: --tree, --table, --metadata, and --output_dir")
}

# Create the base output directory if it doesn't exist
output_dir <- opt$output_dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message(paste("Output directory", output_dir, "created."))
} else {
  message(paste("Output directory", output_dir, "already exists."))
}

# List of subdirectories to create within the output directory
subdirs <- c("checkpoints", "pcoas", "rf", "permanova", "correction")

# Loop through and create subdirectories if they don't exist
for (subdir in subdirs) {
  subdir_path <- file.path(output_dir, subdir)
  if (!dir.exists(subdir_path)) {
    dir.create(subdir_path)
    message(paste("Directory", subdir_path, "created."))
  } else {
    message(paste("Directory", subdir_path, "already exists."))
  }
}

# Assign number of cores, ensuring it is at least 1
num_cores <- max(4, parallel::detectCores() - 10)

# Import phylo tree
tree_qza <- read_qza(opt$tree)
phylo_tree <- tree_qza$data

# Load the feature table in QIIME2 format
raw_taxa_qza <- read_qza(opt$table)

# Extract the feature table data from the QIIME2 object and convert to a data frame
raw_taxa <- as.data.frame(raw_taxa_qza$data)
tr_taxa=t(raw_taxa) # transpose
taxa <- tr_taxa[rowSums(tr_taxa) >= 10000, ] # remove samples below our rarefying depth
taxa <- taxa[, colSums(taxa) != 0] # make sure that features with zero counts are removed

# Load full metadata
metadata=read.table(opt$metadata, header = TRUE, sep = "\t", row.names = 1)

# Find the common row names between taxa and metadata
common_rows <- intersect(rownames(taxa), rownames(metadata))

# Subset both taxa and metadata to keep only rows with common row names
taxa <- taxa[common_rows, ]
metadata <- metadata[common_rows, ]
taxa <- taxa[, colSums(taxa) != 0] # confirm no zero sum features are removed

# Check if both filtered tables have the same number of rows
if (nrow(taxa) == nrow(metadata)) {

} else {
  stop("ERROR: The number of rows in taxa and metadata is not equal")
}

# Assign batch id for correction
batchid = metadata[, c('study_id')]
batchid=factor(batchid)
batchid=droplevels(batchid)


# Assign what studies to use as a reference
diet_counts = tapply(metadata$diet_condensed, metadata$study_id,
                     function(x) length(unique(x)))

batchid_list = names(diet_counts[diet_counts >= 2])
#summary(batchid)

# Assign cavariates of interest
covar = metadata[, c('diet_condensed'), drop = FALSE] # maintains correct structure for the permanova
covar$diet_condensed=as.factor(covar$diet_condensed)
covar=droplevels(covar)
#summary(covar)

save(list = ls(), file = file.path(output_dir, "checkpoints", "preload.RData"))


############################################################################################
#
# Loop through each study as the reference while saving the outputs
#
############################################################################################
cat("\nBeginning 'Each Study As Reference' Analysis\n")
tic("'Each Study As Reference' Analysis Complete")

all_cross_entropy_key_df <- data.frame()
all_cross_entropy_batch_df <- data.frame()

for (study in batchid_list) {

    # Define the checkpoint file path
    checkpoint_file <- file.path(output_dir, "checkpoints", paste0(study, "_loop.RData"))
  
  # Skip this study if the checkpoint already exists
  if (file.exists(checkpoint_file)) {
    cat(paste0("Skipping ", study, " — checkpoint already exists.\n"))
    next  # move on to the next study in the loop
  }


  cat(paste("Beginning Default ConQuR with Ref. Study", study, "\n"))
  tic(paste("Completed Default ConQuR with Ref. Study", study))
  options(warn=-1)
  taxa_corrected = ConQuR(num_core = num_cores,
                          tax_tab=taxa, batchid=batchid,
                          covariates=covar, batch_ref=study)
  assign(paste0("taxa_corrected_ref_", study), taxa_corrected)
  toc()

  cat(paste("Beginning Penalized ConQuR with Ref. Study", study, "\n"))
  tic(paste("Completed Penalized ConQuR with Ref. Study", study))
  options(warn=-1)
  taxa_pen_corrected = ConQuR(num_core = num_cores,
                              tax_tab=taxa, batchid=batchid,
                              covariates=covar, batch_ref=study,
                              logistic_lasso=T, quantile_type="lasso", interplt=T)
  assign(paste0("taxa_pen_corrected_ref_", study), taxa_pen_corrected)
  toc()

  # Saving corrected feature tables
  taxa_corrected_tr = t(taxa_corrected)
  taxa_pen_corrected_tr = t(taxa_pen_corrected)

  # Convert to a data frame and assign the first column name as "SampleID"
  taxa_corrected_tr <- as.data.frame(taxa_corrected_tr)
  taxa_corrected_tr <- cbind(SampleID = rownames(taxa_corrected_tr), taxa_corrected_tr)

  taxa_pen_corrected_tr <- as.data.frame(taxa_pen_corrected_tr)
  taxa_pen_corrected_tr <- cbind(SampleID = rownames(taxa_pen_corrected_tr), taxa_pen_corrected_tr)

  # Reset row names since they are now moved to the first column
  rownames(taxa_corrected_tr) <- NULL
  rownames(taxa_pen_corrected_tr) <- NULL

  # Save table for Qiime2 import
  write.table(taxa_corrected_tr, file = file.path(output_dir, "correction", paste0("taxa_table_corrected_", study, ".tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Default correction feature table saved in correction/\n")

  write.table(taxa_pen_corrected_tr, file = file.path(output_dir, "correction", paste0("taxa_table_penalized_corrected_", study, ".tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Penalized correction feature table saved in correction/\n")

  # Plotting PCoAs
  cat(paste("Beginning Penalized PCoA Plots with Ref. Study", study, "\n"))
  tic(paste("Completed PCoA Plots with Ref. Study", study))

  # Set the output file for the plot
  png(filename = file.path(output_dir, "pcoas", paste0(study, "_pcoa.png")), width = 2400, height = 1600, res = 300)

  par(mfrow=c(3, 3))

  Plot_PCoA(TAX=taxa, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_1", tree=phylo_tree, main=paste("No Correction, WUniFrac: ", study))
  Plot_PCoA(TAX=taxa_corrected, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_1", tree=phylo_tree, main=paste("ConQuR (Default), WUniFrac: ", study))
  Plot_PCoA(TAX=taxa_pen_corrected, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_1", tree=phylo_tree, main=paste("ConQuR (Penalized), WUniFrac: ", study))

  Plot_PCoA(TAX=taxa, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_0.5", tree=phylo_tree, main=paste("No Correction, GUniFrac: ", study))
  Plot_PCoA(TAX=taxa_corrected, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_0.5", tree=phylo_tree, main=paste("ConQuR (Default), GUniFrac: ", study))
  Plot_PCoA(TAX=taxa_pen_corrected, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_0.5", tree=phylo_tree, main=paste("ConQuR (Penalized), GUniFrac: ", study))

  Plot_PCoA(TAX=taxa, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_0", tree=phylo_tree, main=paste("No Correction, UWUniFrac: ", study))
  Plot_PCoA(TAX=taxa_corrected, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_0", tree=phylo_tree, main=paste("ConQuR (Default), UWUniFrac: ", study))
  Plot_PCoA(TAX=taxa_pen_corrected, factor=batchid, dissimilarity="GUniFrac", GUniFrac_type="d_0", tree=phylo_tree, main=paste("ConQuR (Penalized), UWUniFrac: ", study))

  # Close the graphics device
  dev.off()
  toc()

  ### Random forest
  # Whether the corrected taxa read count table can better predict the key variable
  cat(paste("Beginning Random Forest Prediction of Key Variable with Ref. Study", study, "\n"))
  tic(paste("Completed Random Forest Prediction of Key Variable with Ref. Study", study))
  diet_condensed = covar[, 'diet_condensed']
  # Check if 'diet_condensed' is a factor, stop if it's not (has to be a factor to get stratified)
  if (!is.factor(diet_condensed)) {
    stop("'diet_condensed' is not a factor. Please ensure it is a factor to proceed.")
  }
  taxa_result = list(taxa, taxa_corrected, taxa_pen_corrected)

  # Calculate entropy for each taxa correction type
  cross_entropy_results = matrix(ncol=3, nrow=10)
  colnames(cross_entropy_results) = c("No Correction", "ConQuR (Default)", "ConQuR (Penalized)")

  # Set up stratified cross-validation with 5 folds
  set.seed(123)  # Setting seed for reproducibility
  train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = multiClassSummary, savePredictions = TRUE)

  # Loop through each version of the taxa correction types
  for (ii in 1:3) {
    # Prepare the data frame for caret training: taxa_result[ii] and diet_condensed combined
    taxa_data <- data.frame(taxa_result[[ii]])
    taxa_data$diet_condensed <- diet_condensed

    # Train the Random Forest model
    rf_model <- train(
      diet_condensed ~ .,                # Formula: predict diet_condensed using all other columns
      data = taxa_data,                  # Data frame containing features and target
      method = "rf",                     # Random Forest method
      trControl = train_control,         # Cross-validation settings
      importance = TRUE                  # Request feature importance
    )

    # Extract cross-entropy (log loss) from each fold of cross-validation
    fold_results <- rf_model$pred  # Get predictions for all folds
    log_loss_fold <- sapply(split(fold_results, fold_results$Resample), function(fold) {
      predicted_probs <- as.matrix(fold[, levels(diet_condensed)])
      true_labels <- fold$obs
      predicted_probs <- pmax(predicted_probs, 1e-6)
      log_loss <- -mean(log(predicted_probs[cbind(1:nrow(fold), as.numeric(true_labels))]))
      return(log_loss)
    })

    # Store the log loss values (cross-entropy) for each fold
    cross_entropy_results[, ii] <- log_loss_fold
  }
  toc()

  # Convert the cross-entropy results to a data frame
  cat(paste("Beginning Extraction of Random Forest Prediction of Key Variable Results and Plotting with Ref. Study", study, "\n"))
  tic(paste("Completed Extraction of Random Forest Prediction of Key Variable Results and Plotting with Ref. Study", study))
  cross_entropy_key_df <- data.frame(
    Study = rep(study, times = 30),
    Fold = rep(1:10, times = 3),
    Correction = rep(c("No Correction", "ConQuR (Default)", "ConQuR (Penalized)"), each = 10),
    CrossEntropy = as.vector(cross_entropy_results)
  )

  # Set the order of Correction factor levels
  cross_entropy_key_df$Correction <- factor(cross_entropy_key_df$Correction, levels = c("No Correction", "ConQuR (Default)", "ConQuR (Penalized)"))

  # Append cross entropy dataframe to master dataframe
  all_cross_entropy_key_df <- rbind(all_cross_entropy_key_df, cross_entropy_key_df)

  # Perform ANOVA and Tukey
  anova_result <- aov(CrossEntropy ~ Correction, data = cross_entropy_key_df)
  tukey_result <- TukeyHSD(anova_result)
  cld <- multcompLetters4(anova_result, tukey_result)

  # Create the Tk table with factors and 3rd quantile
  Tk <- cross_entropy_key_df %>%
    group_by(Correction) %>%
    dplyr::summarise(mean = mean(CrossEntropy), quant = quantile(CrossEntropy, probs = 0.75)) %>%
    arrange(desc(mean))

  # Extract the compact letter display and add to the Tk table
  cld_df <- as.data.frame.list(cld$Correction)
  Tk$cld <- cld_df$Letters

  # Create the box plot using ggplot2
  ggplot(cross_entropy_key_df, aes(x = Correction, y = CrossEntropy, fill = Correction)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Cross-Log-Loss Variability Across Folds for Different Correction Types",
         x = "Correction Type",
         y = "Cross-Entropy (Log Loss)") +
    theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
    geom_text(data = Tk, aes(x = Correction, y = quant, label = cld), size = 4, vjust=-1, hjust =-1)
  ggsave(file.path(output_dir, "rf", paste0(study, "_key_rf_log_loss_boxplot.png")), width = 8, height = 8, dpi = 300)

  save(list = ls(), file = file.path(output_dir, "checkpoints", paste0(study, "_loop.RData")))

  toc()


}

# Factor order the master dfs
all_cross_entropy_key_df$Correction <- factor(all_cross_entropy_key_df$Correction, levels = c("No Correction", "ConQuR (Default)", "ConQuR (Penalized)"))

save(list = ls(), file = file.path(output_dir, "checkpoints", "all_studies_corrections.RData"))

toc()

load(file.path(output_dir, "checkpoints", "all_studies_corrections.RData"))
############################################################################################
#
# Random Forest Results
#
############################################################################################
cat("Plotting Random Forest Results\n")

### Key prediction
# Filter the all_cross_entropy_key_df table for "No Correction" correction values only
original_values <- all_cross_entropy_key_df %>%
  filter(Correction == "No Correction")

# Split the data by Study to get a list of RF sets for each study
original_values_list <- split(original_values$CrossEntropy, original_values$Study)

# Check if all sets are identical to the first set
all_sets_identical <- all(sapply(original_values_list, function(x) all(x == original_values_list[[1]])))

# Output the result
if (all_sets_identical) {
  cat("All sets of 'No Correction' Entropy values are identical across studies. Proceding to remove duplicates.\n")
  # Remove duplicate "No Correction" rows, keep only one for "No Correction" correction
  rf_combined_unique <- all_cross_entropy_key_df %>%
    filter(!(Correction == "No Correction" & Study != unique(Study)[1])) %>%
    mutate(Study = ifelse(Correction == "No Correction", "No Correction", Study))
} else {
  stop("ERROR: There are differences in the sets of 'No Correction' Entropy values across studies.")
}

# Separate the "No Correction" data from the other corrections
no_correction_data <- rf_combined_unique %>% filter(Correction == "No Correction")
correction_data <- rf_combined_unique %>% filter(Correction != "No Correction")

# Get unique study names excluding "No Correction"
studies_without_nocorr <- unique(correction_data$Study)

# Create color palette
num_studies <- length(unique(batchid)) + 1
cat("There are", num_studies, "studies in this analysis. \n")
palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_studies)

#### Make the above code a variable that counts the number of studies present - have to make your own palette with RColorBrewer


# Plot RF values from all studies
ggplot() +
  # Plot the "No Correction" group with a different box width
  geom_boxplot(data = no_correction_data, aes(x = Correction, y = CrossEntropy, fill = Study), width = 0.33, show.legend = FALSE) +
  # Plot the other groups with standard width
  geom_boxplot(data = correction_data, aes(x = Correction, y = CrossEntropy, fill = Study), width = 1.0) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "solid", color = "black", size = 0.5) +
  theme_classic() +
  labs(
    title = "Log Loss of Predicting Diet Condensed by Correction Type",
    x = "",
    y = "Log Loss",
    fill = "Ref.\nStudy"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(1, "lines")
  ) +
  # Remove "No Correction" from the legend by specifying the exact levels in the breaks
  scale_fill_manual(
    values = palette, 
    breaks = studies_without_nocorr  # This ensures "No Correction" is excluded
  )

# Save the plot to a file
ggsave(file.path(output_dir, "rf", "all_studies_grouped_key_entropy_boxplot.png"), width = 12, height = 8, dpi = 300)

cat("Random Forest Plots Saved to rf/\n")

save(list = ls(), file = file.path(output_dir, "checkpoints", "full_default_program_complete.RData"))

cat("***\nProgram Complete\n***")