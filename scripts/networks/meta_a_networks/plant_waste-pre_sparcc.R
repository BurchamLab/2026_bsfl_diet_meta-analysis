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

dir.create("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/feature_tables/plant_waste/")

set.seed(3)

paste("Working directory has been set.")

############### Load data

metadata = read.delim("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/metadata/bsf_metadata-r.txt", sep = "\t", header = TRUE, row.names = 1)

if ("plant_waste" == "gainesville" || "plant_waste" == "chicken_feed") {
    
    metadata_plant_waste = metadata %>%
  dplyr::filter(diet == "plant_waste")

} else {

metadata_plant_waste = metadata %>%
  dplyr::filter(diet_condensed == "plant_waste")

}

samples_kept = rownames(metadata_plant_waste)

feature_table = read.delim("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/batch_correction/correction/taxa_table_corrected_current_study.tsv", sep = "\t", header = TRUE, row.names = 1)

feature_table_plant_waste = feature_table[, colnames(feature_table) %in% samples_kept]

asvs_kept = rownames(feature_table_plant_waste)

taxa_table = read.table("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/taxonomy/taxa/taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1)

taxa_table_plant_waste = taxa_table[rownames(taxa_table) %in% asvs_kept, ]


### Making a split taxonomy file

sep.table = separate(taxa_table_plant_waste, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", remove = FALSE)

taxonomy_cols = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

for (col in taxonomy_cols) {
  sep.table[[col]] <- sub("^.*__", "", sep.table[[col]])
}

paste("Metadata, feature table, and taxonomy have been loaded.")


# Check dimensions
print(dim(feature_table_plant_waste))  # Should match the number of rows in tax
print(dim(sep.table))    # Ensure tax has the correct number of rows

# Ensure all missing or "NA" strings are set to NA
tax_na_unknown = sep.table
tax_na_unknown[is.na(tax_na_unknown)] <- "Unknown" # Should replace NA values with Unknowns

tax_uncul_unknown = tax_na_unknown
tax_uncul_unknown[tax_uncul_unknown == "uncultured"] <- "Unknown" # Should replace "uncultured" with Unknown

tax = tax_uncul_unknown
tax[tax == ""] = "Unknown"

tax$Taxa_Name = tax$Taxon


# Make color palettes from the full taxa table.

sep.table_full = separate(taxa_table, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", remove = FALSE)

for (col in taxonomy_cols) {
  sep.table_full[[col]] <- sub("^.*__", "", sep.table_full[[col]])
}

sep.table_full[is.na(sep.table_full)] = "Unknown"
sep.table_full[sep.table_full == "uncultured"] = "Unknown"
sep.table_full[sep.table_full == ""] = "Unknown"


unique_genera <- unique(sep.table_full$Genus)
genus_palette <- createPalette(
  length(unique_genera),
  seedcolors = c(
    "#8B0000",  # dark red
    "#B8860B",  # dark goldenrod (yellow/brown)
    "#006400",  # dark green
    "#2F4F4F",  # dark slate gray / teal
    "#00008B"   # dark blue
  ),
    range = c(15, 60)
)
names(genus_palette) <- unique_genera
sep.table_full$Genus_color <- genus_palette[sep.table_full$Genus]
sep.table_full$Genus_color[sep.table_full$Genus == "Unknown"] <- "#555555"

unique_families <- unique(sep.table_full$Family)
family_palette <- createPalette(
  length(unique_families),
  seedcolors = c(
    "#8B0000",  # dark red
    "#B8860B",  # dark goldenrod (yellow/brown)
    "#006400",  # dark green
    "#2F4F4F",  # dark slate gray / teal
    "#00008B"   # dark blue
  ),
  range = c(15, 60)
)
names(family_palette) <- unique_families
sep.table_full$Family_color <- family_palette[sep.table_full$Family]
sep.table_full$Family_color[sep.table_full$Family == "Unknown"] <- "#555555"

unique_orders <- unique(sep.table_full$Order)
order_palette <- createPalette(
  length(unique_genera),
    range = c(15, 60),
    seedcolors = c(
    "#8B0000",  # dark red
    "#B8860B",  # dark goldenrod (yellow/brown)
    "#006400",  # dark green
    "#2F4F4F",  # dark slate gray / teal
    "#00008B"   # dark blue
  )
)
names(order_palette) <- unique_orders
sep.table_full$Order_color <- order_palette[sep.table_full$Order]
sep.table_full$Order_color[sep.table_full$Order == "Unknown"] <- "#555555"

unique_classes <- unique(sep.table_full$Class)
class_palette <- createPalette(
  length(unique_classes),
    range = c(15, 60),
    seedcolors = c(
    "#8B0000",  # dark red
    "#B8860B",  # dark goldenrod (yellow/brown)
    "#006400",  # dark green
    "#2F4F4F",  # dark slate gray / teal
    "#00008B"   # dark blue
  )
)
names(class_palette) <- unique_classes
sep.table_full$Class_color <- class_palette[sep.table_full$Class]
sep.table_full$Class_color[sep.table_full$Class == "Unknown"] <- "#555555"

unique_phyla <- unique(sep.table_full$Phylum)
phylum_palette <- createPalette(
  length(unique_phyla),
    range = c(15, 60),
    seedcolors = c(
    "#8B0000",  # dark red
    "#B8860B",  # dark goldenrod (yellow/brown)
    "#006400",  # dark green
    "#2F4F4F",  # dark slate gray / teal
    "#00008B"   # dark blue
  )
)
names(phylum_palette) <- unique_phyla
sep.table_full$Phylum_color <- phylum_palette[sep.table_full$Phylum]
sep.table_full$Phylum_color[sep.table_full$Phylum == "Unknown"] <- "#555555"


############### Calculate relative abundances

relab.rel <- sweep(feature_table_plant_waste, 2, colSums(feature_table_plant_waste), FUN="/") * 100 # Calculates relative abundances within each sample
 
############### 0.1% abundance filtering

sample_count = nrow(metadata_plant_waste)

mmo <- NULL
mmo$prev <- apply(relab.rel > 0.1, 1, sum, na.rm = T) # Abundance filtering

# Filters data that are not present in at least 10% of samples
relab_abun_filt <- relab.rel[which(mmo$prev >= (sample_count * 0.1)),]
tax_abun_filt <- tax[which(mmo$prev >= (sample_count * 0.1)),]
relab.rel_abun_filt <- relab.rel[which(mmo$prev >= (sample_count * 0.1)),]

paste("Data have been abundance filtered.")

filtered_asvs <- row.names(relab.rel_abun_filt)

filtered_feature_table_plant_waste <- feature_table_plant_waste[filtered_asvs, ]


# Subsample 25 random samples

for (i in 1:100) {

subsamples = sample(colnames(filtered_feature_table_plant_waste), 25)
subsampled_table = filtered_feature_table_plant_waste[, subsamples]

file_path = sprintf("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/feature_tables/plant_waste/plant_waste-table-%d.txt", i)

write.table(x = subsampled_table, file = file_path, sep = "\t", col.names = NA, row.names = TRUE)

}

dir.create("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/plant_waste/")

save.image("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/plant_waste/plant_waste_pre_sparcc.RData")
