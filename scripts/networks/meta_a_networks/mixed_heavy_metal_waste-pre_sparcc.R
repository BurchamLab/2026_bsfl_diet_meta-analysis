# Made by Reese Saho for pre-processing of feature tables before running SparCC on them

invisible(gc()) # clear memory
rm(list = ls()) #clear environment

options(repos = c(CRAN = "https://cloud.r-project.org"))

packages = c("tidyverse", "tools", "igraph", "pals", "reshape2", "tidygraph", "scales", "ggraph", "ggdark", 
             "cowplot", 'dtplyr', 'googledrive', 'googlesheets4', 'haven', 'httr', 'ragg', 'rvest', 'xml2', 'Polychrome', 'GUniFrac', 'compositions', "reticulate")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE) }
  
  library(pkg, character.only = TRUE)
}

paste("Packages have been loaded.")


setwd("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/")


set.seed(3)

paste("Working directory has been set.")

############### Load data

metadata = read.delim("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/metadata/bsf_metadata-r.txt", sep = "\t", header = TRUE, row.names = 1)

metadata_mixed_heavy_metal_waste = metadata %>%
  dplyr::filter(diet_condensed == "mixed_heavy_metal_waste")

samples_kept = rownames(metadata_mixed_heavy_metal_waste)

relab = read.delim("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/batch_correction/correction/taxa_table_corrected_current_study.tsv", sep = "\t", header = TRUE, row.names = 1)

relab_mixed_heavy_metal_waste = relab[, colnames(relab) %in% samples_kept]

asvs_kept = rownames(relab_mixed_heavy_metal_waste)

taxa_table = read.table("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/taxonomy/taxa/taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1)

taxa_table_mixed_heavy_metal_waste = taxa_table[rownames(taxa_table) %in% asvs_kept, ]


### Making a split taxonomy file

sep.table = separate(taxa_table_mixed_heavy_metal_waste, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", remove = FALSE)

taxonomy_cols = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

for (col in taxonomy_cols) {
  #print(sep.table[[col]])
  sep.table[[col]] = substr(sep.table[[col]], 4, nchar(sep.table[[col]]))
}

paste("Metadata, feature table, and taxonomy have been loaded.")


# Check dimensions
print(dim(relab_mixed_heavy_metal_waste))  # Should match the number of rows in tax
print(dim(sep.table))    # Ensure tax has the correct number of rows

# Ensure all missing or "NA" strings are set to NA
tax_na_unknown = sep.table
tax_na_unknown[is.na(tax_na_unknown)] <- "Unknown" # Should replace NA values with Unknowns

tax_uncul_unknown = tax_na_unknown
tax_uncul_unknown[tax_uncul_unknown == "uncultured"] <- "Unknown" # Should replace "uncultured" with Unknown

tax = tax_uncul_unknown

tax$Taxa_Name = tax$Taxon

############### Calculate relative abundances

relab.rel <- sweep(relab_mixed_heavy_metal_waste, 2, colSums(relab_mixed_heavy_metal_waste), FUN="/") * 100 # Calculates relative abundances within each sample
 
############### 0.1% abundance filtering

sample_count = nrow(metadata_mixed_heavy_metal_waste)

mmo <- NULL
mmo$prev <- apply(relab.rel > 0.1, 1, sum, na.rm = T) # Abundance filtering

# Filters data that are not present in at least 20% of samples
relab_abun_filt <- relab_mixed_heavy_metal_waste[which(mmo$prev >= (sample_count * 0.2)),]
tax_abun_filt <- tax[which(mmo$prev >= (sample_count * 0.2)),]
relab.rel_abun_filt <- relab.rel[which(mmo$prev >= (sample_count * 0.2)),]

paste("Data have been abundance filtered.")

# Subsample 25 random samples

for (i in 1:100) {

subsamples = sample(colnames(relab_abun_filt), 25)
subsampled_table = relab_abun_filt[, subsamples]

file_path = sprintf("./mixed_heavy_metal_waste/mixed_heavy_metal_waste-table_%d.txt", i)

write.table(x = subsampled_table, file = file_path, sep = "\t", col.names = NA, row.names = TRUE)

}

save.image("./mixed_heavy_metal_waste/mixed_heavy_metal_waste_pre_sparcc.RData")
