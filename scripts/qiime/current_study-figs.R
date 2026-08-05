# Made by Reese Saho for the purpose of making the alpha diversity, beta diversity, and differential abundance graphs for the Saho BSFL gut samples. 2025/08/03

invisible(gc()) # clear memory
rm(list = ls()) #clear environment

options(repos = c(CRAN = "https://cloud.r-project.org")) # Set a CRAN for package install

# Load in necessary packages

packages = c("devtools",
             "tidyplots",
             "tidyverse",
             "ggplot2", 
             "tidyr",
             "ggpubr",
             "ggrepel")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE) }
  
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

# Set necessary filepaths

# Set the directory
setwd("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa")

### IMPORTING DATA

# Metadata

print("Loading in metadata files.")

q2_metadata = read_q2metadata("./metadata/bsf_metadata-q2r.txt")
metadata = read.table("./metadata/bsf_metadata-r.txt", header = TRUE)

# Feature tables

print("Loading in feature tables.")

pro_table = read.table("./deblur/current_study-pro-table/feature-table-r.tsv", sep = '\t', header = TRUE)
no_pro_table = read.table("./deblur/current_study-no_pro-table/feature-table-r.tsv", sep = '\t', header = TRUE)

# Taxonomies

print("Loading in taxonomy files.")

pro_tax = read.table("./taxonomy/cs-pro-taxa/taxonomy-r.tsv", sep = '\t', header = TRUE)
no_pro_tax = read.table("./taxonomy/cs-no_pro-taxa/taxonomy-r.tsv", sep = '\t', header = TRUE)

# Alpha Diversity Metrics

print("Loading in alpha diversity metrics.")

pro_shannon = read_qza("./core-metrics-phylogeny-current_study-168k-pro/shannon_vector.qza")
no_pro_shannon = read_qza("./core-metrics-phylogeny-current_study-168k-no_pro/shannon_vector.qza")

pro_obs = read_qza("./core-metrics-phylogeny-current_study-168k-pro/observed_features_vector.qza")
no_pro_obs = read_qza("./core-metrics-phylogeny-current_study-168k-no_pro/observed_features_vector.qza")

pro_faith = read_qza("./core-metrics-phylogeny-current_study-168k-pro/faith_pd_vector.qza")
no_pro_faith = read_qza("./core-metrics-phylogeny-current_study-168k-no_pro/faith_pd_vector.qza")

pro_evenness = read_qza("./core-metrics-phylogeny-current_study-168k-pro/evenness_vector.qza")
no_pro_evenness = read_qza("./core-metrics-phylogeny-current_study-168k-no_pro/evenness_vector.qza")

# Beta Diversity Metrics
 
pro_unweighted = read_qza("./core-metrics-phylogeny-current_study-168k-pro/unweighted_unifrac_pcoa_results.qza")
no_pro_unweighted = read_qza("./core-metrics-phylogeny-current_study-168k-no_pro/unweighted_unifrac_pcoa_results.qza")

pro_generalized = read_qza("./core-metrics-phylogeny-current_study-168k-pro/generalized_unifrac_pcoa.qza")
no_pro_generalized = read_qza("./core-metrics-phylogeny-current_study-168k-no_pro/generalized_unifrac_pcoa.qza")

pro_weighted = read_qza("./core-metrics-phylogeny-current_study-168k-pro/weighted_unifrac_pcoa_results.qza")
no_pro_weighted = read_qza("./core-metrics-phylogeny-current_study-168k-no_pro/weighted_unifrac_pcoa_results.qza")

# ANCOM-BC 2 Output Files

pro_lfc = read.csv("./ancom/saho_bsfl-pro-ancom/lfc_slice.csv")
no_pro_lfc = read.csv("./ancom/saho_bsfl-no_pro-ancom/lfc_slice.csv")

pro_q_val = read.csv("./ancom/saho_bsfl-pro-ancom/q_val_slice.csv")
no_pro_q_val = read.csv("./ancom/saho_bsfl-no_pro-ancom/q_val_slice.csv")

# Custom diet labels for the figures

diet_labels = c("tuna" = "Tuna", "bread" = "Bread", "gainesville" = "Gainesville", "spinach" = "Spinach", "coffee_grounds" = "Coffee Grounds", "initial_larvae" = "NA (Initial Larvae)", "cow_manure_probiotic" = "Manure + Probiotics", "cow_manure" = "Manure")

### ALPHA DIVERSITY FIGURES

## Shannon's Richness

# Make the Shannon table match the metadata for merging

pro_shannon = pro_shannon$data %>% rownames_to_column("SampleID")
no_pro_shannon = no_pro_shannon$data %>% rownames_to_column("SampleID")

# Discard samples that aren't present in the Shannon file

pro_shannon_metadata = q2_metadata %>%
  left_join(pro_shannon) %>%
  filter(!is.na(shannon_entropy))
no_pro_shannon_metadata = q2_metadata %>%
  left_join(no_pro_shannon) %>%
  filter(!is.na(shannon_entropy))

# Make the Shannon figure

pro_shannon_metadata %>%
  ggplot(aes(x = diet, y = shannon_entropy)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Shannon Diversity vs. Diet", x = "Diet", y = "Shannon Entropy")

ggsave("./figures/saho_pro_shannon_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_pro_shannon_diet.png", height = 3, width = 4, device = "png")

no_pro_shannon_metadata %>%
  ggplot(aes(x = diet, y = shannon_entropy)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Shannon Diversity vs. Diet", x = "Diet", y = "Shannon Entropy") 
  
ggsave("./figures/saho_no_pro_shannon_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_no_pro_shannon_diet.png", height = 3, width = 4, device = "png")

## Faith's PD

# Make the Faith table match the metadata for merging

pro_faith = pro_faith$data %>% rownames_to_column("SampleID")
no_pro_faith = no_pro_faith$data %>% rownames_to_column("SampleID")

# Discard samples that aren't present in the Faith file

pro_faith_metadata = q2_metadata %>%
  left_join(pro_faith) %>%
  filter(!is.na(faith_pd))
no_pro_faith_metadata = q2_metadata %>%
  left_join(no_pro_faith) %>%
  filter(!is.na(faith_pd))

# Make the Faith figure

pro_faith_metadata %>%
  ggplot(aes(x = diet, y = faith_pd)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Faith's PD vs. Diet", x = "Diet", y = "Faith's PD") 
  
ggsave("./figures/saho_pro_faith_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_pro_faith_diet.png", height = 3, width = 4, device = "png")

no_pro_faith_metadata %>%
  ggplot(aes(x = diet, y = faith_pd)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Faith's PD vs. Diet", x = "Diet", y = "Faith's PD") 

ggsave("./figures/saho_no_pro_faith_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_no_pro_faith_diet.png", height = 3, width = 4, device = "png")

## Observed Features

# Make the Obs table match the metadata for merging

pro_obs = pro_obs$data %>% rownames_to_column("SampleID")
no_pro_obs = no_pro_obs$data %>% rownames_to_column("SampleID")

# Discard samples that aren't present in the Obs file

pro_obs_metadata = q2_metadata %>%
  left_join(pro_obs) %>%
  filter(!is.na(observed_features))
no_pro_obs_metadata = q2_metadata %>%
  left_join(no_pro_obs) %>%
  filter(!is.na(observed_features))

# Make the Obs figure

pro_obs_metadata %>%
  ggplot(aes(x = diet, y = observed_features)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Observed Features vs. Diet", x = "Diet", y = "Observed Features")
  
ggsave("./figures/saho_pro_obs_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_pro_obs_diet.png", height = 3, width = 4, device = "png")

no_pro_obs_metadata %>%
  ggplot(aes(x = diet, y = observed_features)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Observed Features vs. Diet", x = "Diet", y = "Observed Features") 

ggsave("./figures/saho_no_pro_obs_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_no_pro_obs_diet.png", height = 3, width = 4, device = "png")

## Pielou's Evenness

# Make the Evenness table match the metadata for merging

pro_evenness = pro_evenness$data %>% rownames_to_column("SampleID")
no_pro_evenness = no_pro_evenness$data %>% rownames_to_column("SampleID")

# Discard samples that aren't present in the Evenness file

pro_evenness_metadata = q2_metadata %>%
  left_join(pro_evenness) %>%
  filter(!is.na(pielou_evenness))
no_pro_evenness_metadata = q2_metadata %>%
  left_join(no_pro_evenness) %>%
  filter(!is.na(pielou_evenness))

# Make the Evenness figure

pro_evenness_metadata %>%
  ggplot(aes(x = diet, y = pielou_evenness)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Evenness vs. Diet", x = "Diet", y = "Evenness")
  
ggsave("./figures/saho_pro_evenness_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_pro_evenness_diet.png", height = 3, width = 4, device = "png")

no_pro_evenness_metadata %>%
  ggplot(aes(x = diet, y = pielou_evenness)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_x_discrete(labels = diet_labels) +
  theme_q2r() + # There are other themes like theme_bw() or theme_classic()
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5)) +
  labs(title = "Evenness vs. Diet", x = "Diet", y = "Evenness") 

ggsave("./figures/saho_no_pro_evenness_diet.pdf", height = 3, width = 4, device = "pdf") 
ggsave("./figures/saho_no_pro_evenness_diet.png", height = 3, width = 4, device = "png")

### BETA DIVERSITY FIGURES

# Unweighted UniFracs

pro_unweighted$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(q2_metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = diet)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_color_discrete(name="Diet") +
  labs(title = "Unweighted Unifrac With Probiotics") +
  stat_ellipse(aes(group = diet), type = "t", level = 0.95, linetype = 2) +
  coord_fixed()

ggsave("./figures/saho_pro_uw_unifrac.pdf", height=4, width=5, device="pdf")
ggsave("./figures/saho_pro_uw_unifrac.png", height=4, width=5, device="png")

# Unweighted - only GV and manure

pro_unweighted$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(q2_metadata) %>%
  filter(diet %in% c("gainesville", "cow_manure")) %>%
  filter(!SampleID %in% c("InitialLarvae1", "InitialLarvae2", "InitialLarvae3")) %>%
  ggplot(aes(x = PC1, y = PC2, color = diet)) +
  geom_point(alpha=0.5) + 
  stat_ellipse(aes(group = diet), type = "t", level = 0.95, linetype = 2) +
  theme_q2r() +
  scale_color_discrete(name="Diet") +
  labs(title = "Unweighted Unifrac - Gainesville & Manure Samples") +
  coord_fixed()

ggsave("./figures/saho_pro_uw_unifrac_gv_man.pdf", height=4, width=5, device="pdf")
ggsave("./figures/saho_pro_uw_unifrac_gv_man.png", height=4, width=5, device="png")


no_pro_unweighted$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(q2_metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = diet)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_color_discrete(name="Diet") +
  labs(title = "Unweighted Unifrac - No Probiotics") +
  stat_ellipse(aes(group = diet), type = "t", level = 0.95, linetype = 2) +
  coord_fixed()

ggsave("./figures/saho_no_pro_uw_unifrac.pdf", height=4, width=5, device="pdf")
ggsave("./figures/saho_no_pro_uw_unifrac.png", height=4, width=5, device="png")

# Weighted UniFracs

pro_weighted$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(q2_metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = diet)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_color_discrete(name="Diet") +
  labs(title = "Weighted Unifrac With Probiotics") +
  stat_ellipse(aes(group = diet), type = "t", level = 0.95, linetype = 2) +
  coord_fixed()

ggsave("./figures/saho_pro_w_unifrac.pdf", height=4, width=5, device="pdf")
ggsave("./figures/saho_pro_w_unifrac.png", height=4, width=5, device="png")

no_pro_weighted$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(q2_metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = diet)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_color_discrete(name="Diet") +
  labs(title = "Weighted Unifrac - No Probiotics") +
  stat_ellipse(aes(group = diet), type = "t", level = 0.95, linetype = 2) +
  coord_fixed()

ggsave("./figures/saho_no_pro_w_unifrac.pdf", height=4, width=5, device="pdf")
ggsave("./figures/saho_no_pro_w_unifrac.png", height=4, width=5, device="png")
 
# Generalized UniFracs

pro_generalized$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(q2_metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = diet)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_color_discrete(name="Diet") +
  labs(title = "Generalized Unifrac With Probiotics") +
  stat_ellipse(aes(group = diet), type = "t", level = 0.95, linetype = 2) +
  coord_fixed()

ggsave("./figures/saho_pro_g_unifrac.pdf", height=4, width=5, device="pdf")
ggsave("./figures/saho_pro_g_unifrac.png", height=4, width=5, device="png")

no_pro_generalized$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(q2_metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = diet)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_color_discrete(name="Diet") +
  labs(title = "Generalized Unifrac - No Probiotics") +
  stat_ellipse(aes(group = diet), type = "t", level = 0.95, linetype = 2) +
  coord_fixed()

ggsave("./figures/saho_no_pro_g_unifrac.pdf", height=4, width=5, device="pdf")
ggsave("./figures/saho_no_pro_g_unifrac.png", height=4, width=5, device="png")

### DIFFERENTIAL ABUNDANCE

# Making split taxonomy files for labeling

pro_sep_table = separate(pro_tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = FALSE)

taxonomy_cols = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for (col in taxonomy_cols) {
  pro_sep_table[[col]] = substr(pro_sep_table[[col]], 5, nchar(pro_sep_table[[col]]))
}

no_pro_sep_table = separate(no_pro_tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = FALSE)

taxonomy_cols = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for (col in taxonomy_cols) {
  no_pro_sep_table[[col]] = substr(no_pro_sep_table[[col]], 5, nchar(no_pro_sep_table[[col]]))
}

## BREAD VS. GAINESVILLE

bread_table = data.frame(ASV = no_pro_lfc$id, lfc = no_pro_lfc$dietbread, q_val = no_pro_q_val$dietbread)
bread_table$Genus = no_pro_sep_table$Genus[match(bread_table$ASV, no_pro_sep_table$Feature_ID)]
bread_table$Genus[is.na(bread_table$Genus)] = "Unknown"
bread_table$Genus[bread_table$Genus == "uncultured"] = "Unknown"
bread_table$Family = no_pro_sep_table$Family[match(bread_table$ASV, no_pro_sep_table$Feature_ID)]
bread_table$Family[is.na(bread_table$Family)] = "Unknown"
bread_table$log_q_value = -log2(bread_table$q_val)
bread_table$log_q_value[is.infinite(bread_table$log_q_value)] = 985.12322
bread_table$group <- "Not significant"
bread_table$group[bread_table$lfc > 1 & bread_table$log_q_value > 1.3] <- "Inc. DA in Bread"
bread_table$group[bread_table$lfc < -1 & bread_table$log_q_value > 1.3] <- "Inc. DA in Gainesville"
bread_table$Genus[bread_table$Genus == ""] = "Unknown"
bread_table$Family[bread_table$Family == ""] = "Unknown"


# Making the volcano plot - family level labels

ggplot(bread_table, aes(x = lfc, y = log_q_value, , color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Bread" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Bread vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(bread_table, lfc > 2.7 & log_q_value > 128 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(bread_table, lfc < -7.15 & log_q_value > 985.123219 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 20, size = 3, nudge_y = -5, nudge_x = 2
  )

ggsave("./figures/bsf_bread_gv_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_bread_gv_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(bread_table, aes(x = lfc, y = log_q_value, , color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Bread" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Bread vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(bread_table, lfc > 2.7 & log_q_value > 128 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(bread_table, lfc < -7.098 & log_q_value > 985.123219 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_bread_gv_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_bread_gv_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_bread_taxa = bread_table %>%
  filter(q_val < 0.05)

sig_bread_taxa = subset(sig_bread_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_bread = sig_bread_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_bread = sig_bread_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_bread = bind_rows(top_up_bread, top_down_bread)

# Make the bar plot
top_taxa_bread$Genus = make.unique(as.character(top_taxa_bread$Genus))
top_taxa_bread$Genus = factor(top_taxa_bread$Genus, levels = top_taxa_bread$Genus[order(top_taxa_bread$lfc)])

ggplot(top_taxa_bread, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Bread")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_bread_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_bread_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_bread_taxa = bread_table %>%
  filter(q_val < 0.05)

sig_bread_taxa = subset(sig_bread_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_bread = sig_bread_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_bread = sig_bread_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_bread = bind_rows(top_up_bread, top_down_bread)

# Make the bar plot
top_taxa_bread$Family = make.unique(as.character(top_taxa_bread$Family))
top_taxa_bread$Family = factor(top_taxa_bread$Family, levels = top_taxa_bread$Family[order(top_taxa_bread$lfc)])

ggplot(top_taxa_bread, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Bread")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_bread_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_bread_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## COFFEE VS. GAINESVILLE

# Read the ANCOM output files and make a new table that is specific to the diet

coffee_grounds_table = data.frame(ASV = no_pro_lfc$id, lfc = no_pro_lfc$dietcoffee_grounds, q_val = no_pro_q_val$dietcoffee_grounds)
coffee_grounds_table$Genus = no_pro_sep_table$Genus[match(coffee_grounds_table$ASV, no_pro_sep_table$Feature_ID)]
coffee_grounds_table$Genus[is.na(coffee_grounds_table$Genus)] = "Unknown"
coffee_grounds_table$Genus[coffee_grounds_table$Genus == "uncultured"] = "Unknown"
coffee_grounds_table$Family = no_pro_sep_table$Family[match(coffee_grounds_table$ASV, no_pro_sep_table$Feature_ID)]
coffee_grounds_table$Family[is.na(coffee_grounds_table$Family)] = "Unknown"
coffee_grounds_table$Genus[coffee_grounds_table$Genus == ""] = "Unknown"
coffee_grounds_table$Family[coffee_grounds_table$Family == ""] = "Unknown"
coffee_grounds_table$log_q_value = -log2(coffee_grounds_table$q_val)
coffee_grounds_table$log_q_value[is.infinite(coffee_grounds_table$log_q_value)] = 726.0431668
coffee_grounds_table$group <- "Not significant"
coffee_grounds_table$group[coffee_grounds_table$lfc > 1 & coffee_grounds_table$log_q_value > 1.3] <- "Inc. DA in Coffee"
coffee_grounds_table$group[coffee_grounds_table$lfc < -1 & coffee_grounds_table$log_q_value > 1.3] <- "Inc. DA in Gainesville"

# Making the volcano plot - family level labels

ggplot(coffee_grounds_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Coffee" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Coffee vs. Gainesville Differential Abundance - Family-Level") +
  geom_label_repel(
    data = subset(coffee_grounds_table, lfc > 1.5653711 & log_q_value > 8.488717 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(coffee_grounds_table, lfc < -6.90959 & log_q_value > 449.5724853 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_coffee_gv_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_coffee_gv_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(coffee_grounds_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Coffee" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Coffee vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(coffee_grounds_table, lfc > 1.5653711 & log_q_value > 8.488717 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(coffee_grounds_table, lfc < -6.90959 & log_q_value > 449.5724853 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_coffee_gv_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_coffee_gv_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_coffee_grounds_taxa = coffee_grounds_table %>%
  filter(q_val < 0.05)

sig_coffee_grounds_taxa = subset(sig_coffee_grounds_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_coffee_grounds = sig_coffee_grounds_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_coffee_grounds = sig_coffee_grounds_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_coffee_grounds = bind_rows(top_up_coffee_grounds, top_down_coffee_grounds)

# Make the bar plot
top_taxa_coffee_grounds$Genus = make.unique(as.character(top_taxa_coffee_grounds$Genus))
top_taxa_coffee_grounds$Genus = factor(top_taxa_coffee_grounds$Genus, levels = top_taxa_coffee_grounds$Genus[order(top_taxa_coffee_grounds$lfc)])

ggplot(top_taxa_coffee_grounds, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Coffee")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_coffee_grounds_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_coffee_grounds_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_coffee_grounds_taxa = coffee_grounds_table %>%
  filter(q_val < 0.05)

sig_coffee_grounds_taxa = subset(sig_coffee_grounds_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_coffee_grounds = sig_coffee_grounds_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_coffee_grounds = sig_coffee_grounds_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_coffee_grounds = bind_rows(top_up_coffee_grounds, top_down_coffee_grounds)

# Make the bar plot
top_taxa_coffee_grounds$Family = make.unique(as.character(top_taxa_coffee_grounds$Family))
top_taxa_coffee_grounds$Family = factor(top_taxa_coffee_grounds$Family, levels = top_taxa_coffee_grounds$Family[order(top_taxa_coffee_grounds$lfc)])

ggplot(top_taxa_coffee_grounds, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Coffee")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_coffee_grounds_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_coffee_grounds_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## MANURE W/ PROBIOTIC VS. GAINESVILLE

# Read the ANCOM output files and make a new table that is specific to the diet

cow_manure_table = data.frame(ASV = pro_lfc$id, lfc = -(pro_lfc$dietgainesville), q_val = pro_q_val$dietgainesville) # Adding in the -1 because I used the manure as the reference for ANCOM
cow_manure_table$Genus = pro_sep_table$Genus[match(cow_manure_table$ASV, pro_sep_table$Feature_ID)]
cow_manure_table$Genus[is.na(cow_manure_table$Genus)] = "Unknown"
cow_manure_table$Genus[cow_manure_table$Genus == "uncultured"] = "Unknown"
cow_manure_table$Family = pro_sep_table$Family[match(cow_manure_table$ASV, pro_sep_table$Feature_ID)]
cow_manure_table$Family[is.na(cow_manure_table$Family)] = "Unknown"
cow_manure_table$log_q_value = -log2(cow_manure_table$q_val)
cow_manure_table$log_q_value[is.infinite(cow_manure_table$log_q_value)] = 1003.9323
cow_manure_table$group <- "Not significant"
cow_manure_table$group[cow_manure_table$lfc > 1 & cow_manure_table$log_q_value > 1.3] <- "Inc. DA in Probiotic Manure"
cow_manure_table$group[cow_manure_table$lfc < -1 & cow_manure_table$log_q_value > 1.3] <- "Inc. DA in Gainesville"
cow_manure_table$Genus[cow_manure_table$Genus == ""] = "Unknown"
cow_manure_table$Family[cow_manure_table$Family == ""] = "Unknown"

# Making the volcano plot - family level labels

ggplot(cow_manure_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Probiotic Manure" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Manure + Probiotics vs. Gainesville Differential Abundance - Family-Level") +
  geom_label_repel(
    data = subset(cow_manure_table, lfc > 4.48487 & log_q_value > 1003.9322999 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 20, size = 3, nudge_y = -3, nudge_x = 1
  ) +
  geom_label_repel(
    data = subset(cow_manure_table, lfc < -5.7448552 & log_q_value > 481.1879 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -1.2, nudge_x = 0.4
  )

ggsave("./figures/bsf_pro_manure_gv_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_pro_manure_gv_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(cow_manure_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Probiotic Manure" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Manure + Probiotics vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(cow_manure_table, lfc > 4.48487 & log_q_value > 1003.9322999 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 20, size = 3, nudge_y = -5, nudge_x = -3
  ) +
  geom_label_repel(
    data = subset(cow_manure_table, lfc < -5.7448552 & log_q_value > 481.1879 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -5, nudge_x = 0.4
  )

ggsave("./figures/bsf_pro_manure_gv_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_pro_manure_gv_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_cow_manure_taxa = cow_manure_table %>%
  filter(q_val < 0.05)

sig_cow_manure_taxa = subset(sig_cow_manure_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_cow_manure = sig_cow_manure_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_cow_manure = sig_cow_manure_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_cow_manure = bind_rows(top_up_cow_manure, top_down_cow_manure)

# Make the bar plot
top_taxa_cow_manure$Genus = make.unique(as.character(top_taxa_cow_manure$Genus))
top_taxa_cow_manure$Genus = factor(top_taxa_cow_manure$Genus, levels = top_taxa_cow_manure$Genus[order(top_taxa_cow_manure$lfc)])

ggplot(top_taxa_cow_manure, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Probiotic Manure")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_cow_manure_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_cow_manure_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_cow_manure_taxa = cow_manure_table %>%
  filter(q_val < 0.05)

sig_cow_manure_taxa = subset(sig_cow_manure_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_cow_manure = sig_cow_manure_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_cow_manure = sig_cow_manure_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_cow_manure = bind_rows(top_up_cow_manure, top_down_cow_manure)

# Make the bar plot
top_taxa_cow_manure$Family = make.unique(as.character(top_taxa_cow_manure$Family))
top_taxa_cow_manure$Family = factor(top_taxa_cow_manure$Family, levels = top_taxa_cow_manure$Family[order(top_taxa_cow_manure$lfc)])

ggplot(top_taxa_cow_manure, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Probiotic Manure")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_cow_manure_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_cow_manure_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## CONTROL MANURE VS. GAINESVILLE

# Read the ANCOM output files and make a new table that is specific to the diet

manure_control_table = data.frame(ASV = no_pro_lfc$id, lfc = no_pro_lfc$dietcow_manure, q_val = no_pro_q_val$dietcow_manure)
manure_control_table$Genus = no_pro_sep_table$Genus[match(manure_control_table$ASV, no_pro_sep_table$Feature_ID)]
manure_control_table$Genus[is.na(manure_control_table$Genus)] = "Unknown"
manure_control_table$Genus[manure_control_table$Genus == "uncultured"] = "Unknown"
manure_control_table$Family = no_pro_sep_table$Family[match(manure_control_table$ASV, no_pro_sep_table$Feature_ID)]
manure_control_table$Family[is.na(manure_control_table$Family)] = "Unknown"
manure_control_table$log_q_value = -log2(manure_control_table$q_val)
manure_control_table$log_q_value[is.infinite(manure_control_table$log_q_value)] = 952.11724
manure_control_table$group <- "Not significant"
manure_control_table$group[manure_control_table$lfc > 1 & manure_control_table$log_q_value > 1.3] <- "Inc. DA in Manure"
manure_control_table$group[manure_control_table$lfc < -1 & manure_control_table$log_q_value > 1.3] <- "Inc. DA in Gainesville"
manure_control_table$Genus[manure_control_table$Genus == ""] = "Unknown"
manure_control_table$Family[manure_control_table$Family == ""] = "Unknown"

# Making the volcano plot - family level labels

ggplot(manure_control_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Manure" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Manure vs. Gainesville Differential Abundance - Family-Level") +
  geom_label_repel(
    data = subset(manure_control_table, lfc > 3.0201810 & log_q_value > 560 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -2, nudge_x = -2
  ) +
  geom_label_repel(
    data = subset(manure_control_table, lfc < -4 & log_q_value > 375 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_man_con_gv_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_man_con_gv_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(manure_control_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Manure" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Manure vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(manure_control_table, lfc > 3.0201810 & log_q_value > 560 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -1.2, nudge_x = -1
  ) +
  geom_label_repel(
    data = subset(manure_control_table, lfc < -4 & log_q_value > 375 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_man_con_gv_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_man_con_gv_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_manure_control_taxa = manure_control_table %>%
  filter(q_val < 0.05)

sig_manure_control_taxa = subset(sig_manure_control_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_manure_control = sig_manure_control_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_manure_control = sig_manure_control_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_manure_control = bind_rows(top_up_manure_control, top_down_manure_control)

# Make the bar plot
top_taxa_manure_control$Genus = make.unique(as.character(top_taxa_manure_control$Genus))
top_taxa_manure_control$Genus = factor(top_taxa_manure_control$Genus, levels = top_taxa_manure_control$Genus[order(top_taxa_manure_control$lfc)])

ggplot(top_taxa_manure_control, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Manure")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_manure_control_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_manure_control_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_manure_control_taxa = manure_control_table %>%
  filter(q_val < 0.05)

sig_manure_control_taxa = subset(sig_manure_control_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_manure_control = sig_manure_control_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_manure_control = sig_manure_control_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_manure_control = bind_rows(top_up_manure_control, top_down_manure_control)

# Make the bar plot
top_taxa_manure_control$Family = make.unique(as.character(top_taxa_manure_control$Family))
top_taxa_manure_control$Family = factor(top_taxa_manure_control$Family, levels = top_taxa_manure_control$Family[order(top_taxa_manure_control$lfc)])

ggplot(top_taxa_manure_control, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Manure")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_manure_control_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_manure_control_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

#### SPINACH VS. GAINESVILLE

# Read the ANCOM output files and make a new table that is specific to the diet

spinach_table = data.frame(ASV = no_pro_lfc$id, lfc = no_pro_lfc$dietspinach, q_val = no_pro_q_val$dietspinach)
spinach_table$Genus = no_pro_sep_table$Genus[match(spinach_table$ASV, no_pro_sep_table$Feature_ID)]
spinach_table$Genus[is.na(spinach_table$Genus)] = "Unknown"
spinach_table$Genus[spinach_table$Genus == "uncultured"] = "Unknown"
spinach_table$Family = no_pro_sep_table$Family[match(spinach_table$ASV, no_pro_sep_table$Feature_ID)]
spinach_table$Family[is.na(spinach_table$Family)] = "Unknown"
spinach_table$log_q_value = -log2(spinach_table$q_val)
spinach_table$log_q_value[is.infinite(spinach_table$log_q_value)] = 959.92040
spinach_table$group <- "Not significant"
spinach_table$group[spinach_table$lfc > 1 & spinach_table$log_q_value > 1.3] <- "Inc. DA in Spinach"
spinach_table$group[spinach_table$lfc < -1 & spinach_table$log_q_value > 1.3] <- "Inc. DA in Gainesville"
spinach_table$Genus[spinach_table$Genus == ""] = "Unknown"
spinach_table$Family[spinach_table$Family == ""] = "Unknown"

# Making the volcano plot - family level labels

ggplot(spinach_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Spinach" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Spinach vs. Gainesville Differential Abundance - Family-Level") +
  geom_label_repel(
    data = subset(spinach_table, lfc > 2.26400 & log_q_value > 65.7891 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(spinach_table, lfc < -4.949252 & log_q_value > 630.7659 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_spinach_gv_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_spinach_gv_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(spinach_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Spinach" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Spinach vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(spinach_table, lfc > 2.26400 & log_q_value > 65.7891 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(spinach_table, lfc < -4.949252 & log_q_value > 630.7659 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_spinach_gv_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_spinach_gv_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_spinach_taxa = spinach_table %>%
  filter(q_val < 0.05)

sig_spinach_taxa = subset(sig_spinach_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_spinach = sig_spinach_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_spinach = sig_spinach_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_spinach = bind_rows(top_up_spinach, top_down_spinach)

# Make the bar plot
top_taxa_spinach$Genus = make.unique(as.character(top_taxa_spinach$Genus))
top_taxa_spinach$Genus = factor(top_taxa_spinach$Genus, levels = top_taxa_spinach$Genus[order(top_taxa_spinach$lfc)])

ggplot(top_taxa_spinach, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Spinach")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_spinach_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_spinach_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_spinach_taxa = spinach_table %>%
  filter(q_val < 0.05)

sig_spinach_taxa = subset(sig_spinach_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_spinach = sig_spinach_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_spinach = sig_spinach_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_spinach = bind_rows(top_up_spinach, top_down_spinach)

# Make the bar plot
top_taxa_spinach$Family = make.unique(as.character(top_taxa_spinach$Family))
top_taxa_spinach$Family = factor(top_taxa_spinach$Family, levels = top_taxa_spinach$Family[order(top_taxa_spinach$lfc)])

ggplot(top_taxa_spinach, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Spinach")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_spinach_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_spinach_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

#### TUNA VS. GAINESVILLE

# Read the ANCOM output files and make a new table that is specific to the diet

tuna_table = data.frame(ASV = no_pro_lfc$id, lfc = no_pro_lfc$diettuna, q_val = no_pro_q_val$diettuna)
tuna_table$Genus = no_pro_sep_table$Genus[match(tuna_table$ASV, no_pro_sep_table$Feature_ID)]
tuna_table$Genus[is.na(tuna_table$Genus)] = "Unknown"
tuna_table$Genus[tuna_table$Genus == "uncultured"] = "Unknown"
tuna_table$Family = no_pro_sep_table$Family[match(tuna_table$ASV, no_pro_sep_table$Feature_ID)]
tuna_table$Family[is.na(tuna_table$Family)] = "Unknown"
tuna_table$log_q_value = -log2(tuna_table$q_val)
tuna_table$log_q_value[is.infinite(tuna_table$log_q_value)] = 958.01553
tuna_table$group <- "Not significant"
tuna_table$group[tuna_table$lfc > 1 & tuna_table$log_q_value > 1.3] <- "Inc. DA in Tuna"
tuna_table$group[tuna_table$lfc < -1 & tuna_table$log_q_value > 1.3] <- "Inc. DA in Gainesville"
tuna_table$Genus[tuna_table$Genus == ""] = "Unknown"
tuna_table$Family[tuna_table$Family == ""] = "Unknown"

# Making the volcano plot - family level labels

ggplot(tuna_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Tuna" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Tuna vs. Gainesville Differential Abundance - Family-Level") +
  geom_label_repel(
    data = subset(tuna_table, lfc > 3.32988 & log_q_value > 958.0155 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(tuna_table, lfc < -6.672 & log_q_value > 958.0155 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_tuna_gv_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_tuna_gv_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(tuna_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Tuna" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Tuna vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(tuna_table, lfc > 3.32988 & log_q_value > 958.0155 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(tuna_table, lfc < -6.672 & log_q_value > 958.0155 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 0.6, nudge_x = 0.4
  )

ggsave("./figures/bsf_tuna_gv_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_tuna_gv_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_tuna_taxa = tuna_table %>%
  filter(q_val < 0.05)

sig_tuna_taxa = subset(sig_tuna_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_tuna = sig_tuna_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_tuna = sig_tuna_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_tuna = bind_rows(top_up_tuna, top_down_tuna)

# Make the bar plot
top_taxa_tuna$Genus = make.unique(as.character(top_taxa_tuna$Genus))
top_taxa_tuna$Genus = factor(top_taxa_tuna$Genus, levels = top_taxa_tuna$Genus[order(top_taxa_tuna$lfc)])

ggplot(top_taxa_tuna, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Tuna")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_tuna_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_tuna_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_tuna_taxa = tuna_table %>%
  filter(q_val < 0.05)

sig_tuna_taxa = subset(sig_tuna_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_tuna = sig_tuna_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_tuna = sig_tuna_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_tuna = bind_rows(top_up_tuna, top_down_tuna)

# Make the bar plot
top_taxa_tuna$Family = make.unique(as.character(top_taxa_tuna$Family))
top_taxa_tuna$Family = factor(top_taxa_tuna$Family, levels = top_taxa_tuna$Family[order(top_taxa_tuna$lfc)])

ggplot(top_taxa_tuna, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Tuna")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_tuna_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_tuna_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

save.image("saho1.RData")


#### INTIAL LARVAE VS. GAINESVILLE

# Read the ANCOM output files and make a new table that is specific to the diet

initial_larvae_table = data.frame(ASV = no_pro_lfc$id, lfc = no_pro_lfc$dietinitial_larvae, q_val = no_pro_q_val$dietinitial_larvae)
initial_larvae_table$Genus = no_pro_sep_table$Genus[match(initial_larvae_table$ASV, no_pro_sep_table$Feature_ID)]
initial_larvae_table$Genus[is.na(initial_larvae_table$Genus)] = "Unknown"
initial_larvae_table$Genus[initial_larvae_table$Genus == "uncultured"] = "Unknown"
initial_larvae_table$Family = no_pro_sep_table$Family[match(initial_larvae_table$ASV, no_pro_sep_table$Feature_ID)]
initial_larvae_table$Family[is.na(initial_larvae_table$Family)] = "Unknown"
initial_larvae_table$log_q_value = -log2(initial_larvae_table$q_val)
initial_larvae_table$log_q_value[is.infinite(initial_larvae_table$log_q_value)] = 904.98975
initial_larvae_table$group <- "Not significant"
initial_larvae_table$group[initial_larvae_table$lfc > 1 & initial_larvae_table$log_q_value > 1.3] <- "Inc. DA in Initial larvae"
initial_larvae_table$group[initial_larvae_table$lfc < -1 & initial_larvae_table$log_q_value > 1.3] <- "Inc. DA in Gainesville"
initial_larvae_table$Genus[initial_larvae_table$Genus == ""] = "Unknown"
initial_larvae_table$Family[initial_larvae_table$Family == ""] = "Unknown"

# Making the volcano plot - family level labels

ggplot(initial_larvae_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Initial larvae" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Initial Larvae vs. Gainesville Differential Abundance - Family-Level") +
  geom_label_repel(
    data = subset(initial_larvae_table, lfc > 5.409486 & log_q_value > 904.9897 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -10, nudge_x = -3
  ) +
  geom_label_repel(
    data = subset(initial_larvae_table, lfc < -6.672 & log_q_value > 904.9897 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -3, nudge_x = 0.4
  )

ggsave("./figures/bsf_initial_larvae_gv_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_initial_larvae_gv_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(initial_larvae_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Initial larvae" = "red",
      "Inc. DA in Gainesville" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Initial Larvae vs. Gainesville Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(initial_larvae_table, lfc > 5.409486 & log_q_value > 904.9897 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -50, nudge_x = -1
  ) +
  geom_label_repel(
    data = subset(initial_larvae_table, lfc < -6.672 & log_q_value > 904.9897 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -50, nudge_x = 4
  )

ggsave("./figures/bsf_initial_larvae_gv_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_initial_larvae_gv_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_initial_larvae_taxa = initial_larvae_table %>%
  filter(q_val < 0.05)

sig_initial_larvae_taxa = subset(sig_initial_larvae_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_initial_larvae = sig_initial_larvae_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_initial_larvae = sig_initial_larvae_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_initial_larvae = bind_rows(top_up_initial_larvae, top_down_initial_larvae)

# Make the bar plot
top_taxa_initial_larvae$Genus = make.unique(as.character(top_taxa_initial_larvae$Genus))
top_taxa_initial_larvae$Genus = factor(top_taxa_initial_larvae$Genus, levels = top_taxa_initial_larvae$Genus[order(top_taxa_initial_larvae$lfc)])

ggplot(top_taxa_initial_larvae, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Initial Larvae")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_initial_larvae_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_initial_larvae_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_initial_larvae_taxa = initial_larvae_table %>%
  filter(q_val < 0.05)

sig_initial_larvae_taxa = subset(sig_initial_larvae_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_initial_larvae = sig_initial_larvae_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_initial_larvae = sig_initial_larvae_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_initial_larvae = bind_rows(top_up_initial_larvae, top_down_initial_larvae)

# Make the bar plot
top_taxa_initial_larvae$Family = make.unique(as.character(top_taxa_initial_larvae$Family))
top_taxa_initial_larvae$Family = factor(top_taxa_initial_larvae$Family, levels = top_taxa_initial_larvae$Family[order(top_taxa_initial_larvae$lfc)])

ggplot(top_taxa_initial_larvae, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Gainesville", "Inc. DA in Initial Larvae")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_initial_larvae_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_initial_larvae_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## MANURE W/ PROBIOTIC VS. CONTROL MANURE

# Read the ANCOM output files and make a new table that is specific to the diet

manure_comp_table = data.frame(ASV = pro_lfc$id, lfc = -(pro_lfc$dietcow_manure), q_val = pro_q_val$dietcow_manure) # Adding in the -1 because I used the probiotic as the reference for ANCOM
manure_comp_table$Genus = pro_sep_table$Genus[match(manure_comp_table$ASV, pro_sep_table$Feature_ID)]
manure_comp_table$Genus[is.na(manure_comp_table$Genus)] = "Unknown"
manure_comp_table$Genus[manure_comp_table$Genus == "uncultured"] = "Unknown"
manure_comp_table$Family = pro_sep_table$Family[match(manure_comp_table$ASV, pro_sep_table$Feature_ID)]
manure_comp_table$Family[is.na(manure_comp_table$Family)] = "Unknown"
manure_comp_table$log_q_value = -log2(manure_comp_table$q_val)
manure_comp_table$log_q_value[is.infinite(manure_comp_table$log_q_value)] = 636.0380919
manure_comp_table$group <- "Not significant"
manure_comp_table$group[manure_comp_table$lfc > 1 & manure_comp_table$log_q_value > 1.3] <- "Inc. DA in Probiotic Manure"
manure_comp_table$group[manure_comp_table$lfc < -1 & manure_comp_table$log_q_value > 1.3] <- "Inc. DA in Control Manure"
manure_comp_table$Genus[manure_comp_table$Genus == ""] = "Unknown"
manure_comp_table$Family[manure_comp_table$Family == ""] = "Unknown"

# Making the volcano plot - family level labels

ggplot(manure_comp_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Probiotic Manure" = "red",
      "Inc. DA in Control Manure" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Manure + Probiotics vs. Control Manure Differential Abundance - Family-Level") +
  geom_label_repel(
    data = subset(manure_comp_table, lfc > 1.920668 & log_q_value > 61.769198 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 3, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(manure_comp_table, lfc < -1.807680 & log_q_value > 50.815300 & !grepl("Unknown", Genus)),
    aes(label = Family),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 3, nudge_x = 0.4
  )

ggsave("./figures/bsf_manure_comp_diff_abun_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_manure_comp_diff_abun_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

# Making the volcano plot - genus level labels

ggplot(manure_comp_table, aes(x = lfc, y = log_q_value, color = group)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[2]~italic(Q)~adjusted)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    name = "Differential Abundance",
    values = c(
      "Inc. DA in Probiotic Manure" = "red",
      "Inc. DA in Control Manure" = "blue",
      "Not significant" = "lightgrey"
    )
  ) +
  labs( title = "Manure + Probiotics vs. Control Manure Differential Abundance - Genus-Level") +
  geom_label_repel(
    data = subset(manure_comp_table, lfc > 2 & log_q_value > 60.4 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = -1, nudge_x = 0.4
  ) +
  geom_label_repel(
    data = subset(manure_comp_table, lfc < -1.8 & log_q_value > 50 & !grepl("Unknown", Genus)),
    aes(label = Genus),
    fill = "white", color = "black", segment.color = "black",
    box.padding = 0.5, force = 7, size = 3, nudge_y = 3, nudge_x = 0.4
  )

ggsave("./figures/bsf_manure_comp_diff_abun_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_manure_comp_diff_abun_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_manure_comp_taxa = manure_comp_table %>%
  filter(q_val < 0.05)

sig_manure_comp_taxa = subset(sig_manure_comp_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_manure_comp = sig_manure_comp_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_manure_comp = sig_manure_comp_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_manure_comp = bind_rows(top_up_manure_comp, top_down_manure_comp)

# Make the bar plot
top_taxa_manure_comp$Genus = make.unique(as.character(top_taxa_manure_comp$Genus))
top_taxa_manure_comp$Genus = factor(top_taxa_manure_comp$Genus, levels = top_taxa_manure_comp$Genus[order(top_taxa_manure_comp$lfc)])

ggplot(top_taxa_manure_comp, aes(x = Genus, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Control Manure", "Inc. DA in Probiotic Manure")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_manure_comp_gv_barplot_genus_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_manure_comp_gv_barplot_genus_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

## Make a barplot showing the most up- and down-regulated significant taxa in each diet

# Filter out insignificant taxa
sig_manure_comp_taxa = manure_comp_table %>%
  filter(q_val < 0.05)

sig_manure_comp_taxa = subset(sig_manure_comp_taxa, !grepl("Unknown", Genus))

# Grab the most and least differentially abundant taxa
top_up_manure_comp = sig_manure_comp_taxa %>%
  arrange(desc(lfc)) %>%
  head(10)

top_down_manure_comp = sig_manure_comp_taxa %>%
  arrange(lfc) %>%
  head(10)

top_taxa_manure_comp = bind_rows(top_up_manure_comp, top_down_manure_comp)

# Make the bar plot
top_taxa_manure_comp$Family = make.unique(as.character(top_taxa_manure_comp$Family))
top_taxa_manure_comp$Family = factor(top_taxa_manure_comp$Family, levels = top_taxa_manure_comp$Family[order(top_taxa_manure_comp$lfc)])

ggplot(top_taxa_manure_comp, aes(x = Family, y = lfc, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Inc. DA in Control Manure", "Inc. DA in Probiotic Manure")) +
  labs(
    x = NULL,
    y = "Log2 Fold Change",
    title = "Top Differentially Abundant Taxa",
    fill = "Direction"
  ) +
  coord_flip() +
  theme_minimal()

ggsave("./figures/bsf_manure_comp_gv_barplot_family_250803.png", device = png, width = 8, height = 8, unit = "in")
ggsave("./figures/bsf_manure_comp_gv_barplot_family_250803.pdf", device = pdf, width = 8, height = 8, unit = "in")

save.image("./figures/saho_figures.RData")

### SEQUENCING COUNT BOX PLOTS

# Transpose the feature table to calculate sample counts and later match the metadata

t_pro_table = t(pro_table)
t_pro_table = as.data.frame(t_pro_table)
colnames(t_pro_table) = t_pro_table[1,]
t_pro_table = t_pro_table[-1,]
t_pro_table[] <- lapply(t_pro_table, function(x) as.numeric(as.character(x)))

t_pro_table$read_count = rowSums(t_pro_table[, sapply(t_pro_table, is.numeric)])

t_pro_table$SampleID <- rownames(t_pro_table)
feature_table_metadata = left_join(t_pro_table, metadata, by = "SampleID")

# Make the boxplot

ggplot(feature_table_metadata, aes(x = diet, y = read_count)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(aes(y = read_count)) +
  xlab("Diet") +
  ylab("# of Sequencing Reads")

ggsave("./figures/saho_sample_read_count.pdf", device = pdf, width = 8, height = 8, unit = "in")
ggsave("./figures/saho_sample_read_count.png", device = png, width = 8, height = 8, unit = "in")

# Gainesville & Manure relative abundance plots

gv_man_table = pro_table %>%
    filter((diet %in% c("gainesville", "cow_manure"))
    
gv_man <- gv_man_table
gv_man[,-1] <- gv_man[,-1] / rowSums(gv_man[,-1])

# Reshape to long format for ggplot
gv_man_ <- melt(abundance_rel, id.vars = "Sample",
                       variable.name = "Taxon",
                       value.name = "RelativeAbundance")

# Plot
ggplot(abundance_long, aes(x = Sample, y = RelativeAbundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Relative Abundance") +
  theme_minimal()


