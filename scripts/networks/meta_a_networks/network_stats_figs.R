# Made by Reese Saho for making figures for the network stats analyses. 1/26/2026

invisible(gc()) # clear memory
rm(list = ls()) #clear environment


packages = c('dplyr', 'ggplot2', 'tidyr')

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE) }
  
  library(pkg, character.only = TRUE)
}

# Load in the bootstrapped AUC files for each diet

animal_product_metrics = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_product/random_auc_table.tsv", sep = "\t", header = TRUE)
animal_waste_metrics = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/random_auc_table.tsv", sep = "\t", header = TRUE)
mixed_product_metrics = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/mixed_product/random_auc_table.tsv", sep = "\t", header = TRUE)
food_waste_metrics = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/food_waste/random_auc_table.tsv", sep = "\t", header = TRUE)
plant_product_metrics = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/plant_product/random_auc_table.tsv", sep = "\t", header = TRUE)
plant_waste_metrics = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/plant_waste/random_auc_table.tsv", sep = "\t", header = TRUE)

# Load in the global metric table for each diet

animal_product_global = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_product/metric_auc.tsv", sep = "\t", header = TRUE)
animal_waste_global = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/metric_auc.tsv", sep = "\t", header = TRUE)
mixed_product_global = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/mixed_product/metric_auc.tsv", sep = "\t", header = TRUE)
food_waste_global = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/food_waste/metric_auc.tsv", sep = "\t", header = TRUE)
plant_product_global = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/plant_product/metric_auc.tsv", sep = "\t", header = TRUE)
plant_waste_global = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/plant_waste/metric_auc.tsv", sep = "\t", header = TRUE)

setwd("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/")

# Make two data frames with the diet and AUC values

mod_auc_values = data.frame(
  animal_product = animal_product_metrics$AUC_Modularity,
  animal_waste   = animal_waste_metrics$AUC_Modularity,
  mixed_product  = mixed_product_metrics$AUC_Modularity,
  food_waste     = food_waste_metrics$AUC_Modularity,
  plant_product  = plant_product_metrics$AUC_Modularity,
  plant_waste    = plant_waste_metrics$AUC_Modularity
)

clus_auc_values = data.frame(
  animal_product = animal_product_metrics$AUC_Clustering,
  animal_waste   = animal_waste_metrics$AUC_Clustering,
  mixed_product  = mixed_product_metrics$AUC_Clustering,
  mixed_waste    = food_waste_metrics$AUC_Clustering,
  plant_product  = plant_product_metrics$AUC_Clustering,
  plant_waste    = plant_waste_metrics$AUC_Clustering
)

path_len_auc_values = data.frame(
  animal_product = animal_product_metrics$AUC_AvgPathLength,
  animal_waste   = animal_waste_metrics$AUC_AvgPathLength,
  mixed_product  = mixed_product_metrics$AUC_AvgPathLength,
  mixed_waste    = food_waste_metrics$AUC_AvgPathLength,
  plant_product  = plant_product_metrics$AUC_AvgPathLength,
  plant_waste    = plant_waste_metrics$AUC_AvgPathLength
)

# Subtract the original global value from the random removal value

mod_auc_values_norm = data.frame(
  animal_product = (((animal_product_global$Modularity - animal_product_metrics$AUC_Modularity) / animal_product_metrics$AUC_Modularity) * 100),
  animal_waste   = (((animal_waste_global$Modularity - animal_waste_metrics$AUC_Modularity) / animal_waste_metrics$AUC_Modularity) * 100),
  mixed_product  = (((mixed_product_global$Modularity - mixed_product_metrics$AUC_Modularity) / mixed_product_metrics$AUC_Modularity) * 100),
  food_waste     = (((food_waste_global$Modularity - food_waste_metrics$AUC_Modularity) / food_waste_metrics$AUC_Modularity) * 100),
  plant_product  = (((plant_product_global$Modularity - plant_product_metrics$AUC_Modularity) / plant_product_metrics$AUC_Modularity) * 100),
  plant_waste    = (((plant_waste_global$Modularity - plant_waste_metrics$AUC_Modularity) / plant_waste_metrics$AUC_Modularity) * 100)
)

clus_auc_values = data.frame(
  animal_product = animal_product_metrics$AUC_Clustering,
  animal_waste   = animal_waste_metrics$AUC_Clustering,
  mixed_product  = mixed_product_metrics$AUC_Clustering,
  mixed_waste    = food_waste_metrics$AUC_Clustering,
  plant_product  = plant_product_metrics$AUC_Clustering,
  plant_waste    = plant_waste_metrics$AUC_Clustering
)

path_len_auc_values = data.frame(
  animal_product = animal_product_metrics$AUC_AvgPathLength,
  animal_waste   = animal_waste_metrics$AUC_AvgPathLength,
  mixed_product  = mixed_product_metrics$AUC_AvgPathLength,
  mixed_waste    = food_waste_metrics$AUC_AvgPathLength,
  plant_product  = plant_product_metrics$AUC_AvgPathLength,
  plant_waste    = plant_waste_metrics$AUC_AvgPathLength
)

# Convert the tables into ggplot-friendly forms

mod_table = mod_auc_values_norm %>%
  pivot_longer(cols = all_of(names(mod_auc_values)), names_to = "Diet", values_to = "AUC_values")

apl_table = path_len_auc_values %>%
  pivot_longer(cols = all_of(names(path_len_auc_values)), names_to = "Diet", values_to = "AUC_values")

clus_table = clus_auc_values %>%
  pivot_longer(cols = all_of(names(clus_auc_values)), names_to = "Diet", values_to = "AUC_values")

# Do stats testing

kruskal.test(AUC_values ~ Diet, data = mod_table)
kruskal.test(AUC_values ~ Diet, data = clus_table)
kruskal.test(AUC_values ~ Diet, data = apl_table)

pairwise.wilcox.test(
  mod_table$AUC_values,
  mod_table$Diet,
  p.adjust.method = "BH"  # Benjamini-Hochberg FDR correction
)

pairwise.wilcox.test(
  clus_table$AUC_values,
  clus_table$Diet,
  p.adjust.method = "BH"  # Benjamini-Hochberg FDR correction
)

pairwise.wilcox.test(
  apl_table$AUC_values,
  apl_table$Diet,
  p.adjust.method = "BH"  # Benjamini-Hochberg FDR correction
)


# Make box and whisker plots of the distribution of the modularity values

mod_plot = ggplot(mod_table, aes(x = Diet, y = AUC_values)) +
  geom_violin() +
  geom_boxplot(color = "black", alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_jitter(aes(y = AUC_values)) +
  ggtitle("Modularity AUC by Diet") +
  xlab("Diet") +
  ylab("AUC % Change")

ggsave(filename = "mod_boxplot.pdf", plot = mod_plot, device = "pdf", width = 10, height = 8, units = "in")

# Make box and whisker plots of the distribution of the clustering coefficient values

clus_plot = ggplot(clus_table, aes(x = Diet, y = AUC_values)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(aes(y = AUC_values)) +
  ggtitle("Clustering Coefficient AUC by Diet") +
  xlab("Diet") +
  ylab("AUC - Clustering Coefficient")

ggsave(filename = "clus_boxplot.pdf", plot = clus_plot, device = "pdf", width = 10, height = 8, units = "in")

# Make box and whisker plots of the distribution of the path length values

apl_plot = ggplot(apl_table, aes(x = Diet, y = AUC_values)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(aes(y = AUC_values)) +
  ggtitle("Average Path Length AUC by Diet") +
  xlab("Diet") +
  ylab("AUC - Average Path Length")

ggsave(filename = "path_length_boxplot.pdf", plot = apl_plot, device = "pdf", width = 10, height = 8, units = "in")

getwd()
