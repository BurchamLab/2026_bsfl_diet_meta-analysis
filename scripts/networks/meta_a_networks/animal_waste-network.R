# Made by Reese Saho for the purpose of creating networks and edge files for use in figure generation

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

message("Packages have been loaded.")

# Load in the pre-SparCC file

load("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/animal_waste_pre_sparcc.RData")


# Load in all correlation and p-value matrices together as lists

base_path = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/sparcc_output/animal_waste"

cor_files = list.files(path = base_path, pattern = ".*cor_sparcc\\.csv$", 
                       recursive = TRUE, 
                       full.names = TRUE)

cor_data_list <- lapply(cor_files, read.csv, row.names = 1)

message("List of correlation files has been made.")

names(cor_data_list) <- basename(dirname(cor_files))

p_val_files = list.files(path = base_path, pattern = ".*pvals_one_sided\\.csv$", 
                         recursive = TRUE, 
                         full.names = TRUE)

p_val_data_list <- lapply(p_val_files, read.csv, row.names = 1)

names(p_val_data_list) <- basename(dirname(cor_files))


message("List of p-value files has been made.")


# Set up a for loop to run through all of the correlation and p-value matrices

for (name in names(cor_data_list)) {
  cor_mat = cor_data_list[[name]]
  p_val = p_val_data_list[[name]]
  
  p_val_mat = as.matrix(p_val)
  
  cor_mat_mat = as.matrix(cor_mat)
  
  # Remove all correlations that are not statistically significant
  
  cor_mat_mat[p_val_mat >= 0.05] = 0.00001
  
  cor_mat_mat <- ifelse(p_val_mat >= 0.05, 0.00001, cor_mat_mat)
  
  message("Statistically insignificant correlations have been removed.")
  
  # Rename the correlation matrix with the corresponding ASVs
  
  rownames(cor_mat_mat) <- rownames(subsampled_table)
  colnames(cor_mat_mat) <- rownames(subsampled_table)
  
  ############### Link ASVs to taxonomies
  
  asv_seqs = row.names(tax_abun_filt)
  
  ############### Make network
  
  cor.thr <- 0.3
  
  otu_net <- graph_from_adjacency_matrix(cor_mat_mat, weighted = TRUE, mode = "upper", diag = FALSE)
  otu_net <- delete_edges(otu_net, which(E(otu_net)$weight < cor.thr | E(otu_net)$weight %in% NA))
  otu_net <- igraph::simplify(otu_net, remove.multiple = T, remove.loops = T)
  
  message("Network has been made.")
  
  ############### Get Family-rank taxonomy for each node
  
  tax.sel <- tax %>% select("Taxa_Name")
  tax.sel <- tibble::rownames_to_column(tax.sel, "name")
  
  nodes <- data.frame(name = V(otu_net)$name)
  tax.flt <- tax.sel %>% filter(name %in% nodes$name)
  means <- data.frame(name = rownames(relab.rel),
                      meanAbundance = rowMeans(relab.rel))
  tax.mrg <- merge(means, tax.flt, by = "name", sort = FALSE)
  
  message("Taxnonomies have been rank-ordered.")
  
  ###############  Create list of top 15 most abundant taxonomies
  
  tax.top <- tax.mrg %>%
    group_by(.data[["Taxa_Name"]]) %>%
    summarise(sum = sum(meanAbundance)) %>%
    pull("Taxa_Name")
  
  
  ############### Change taxonomy to "Other" for entries outside of the top 15
  final.tax <- case_when(tax.mrg[["Taxa_Name"]] %in% tax.top ~ tax.mrg[["Taxa_Name"]], TRUE ~ "Other")
  
  
  message("Low abundance taxa have been changed to Other.")
  
  
  ############## Rename nodes to numbers
  
  asv_names <- rownames(subsampled_table)
  
  # Map ASVs to numbers
  
  # full_taxonomy = read.table(file = "~/Library/CloudStorage/OneDrive-UniversityofTennessee/UTK/Research/Code/bsf_metaa/SparCC-master/plant/bsf_plant_taxonomy.txt", sep = "\t", header = TRUE)
  
  full_taxonomy = read.table(file = "/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/taxonomy/taxa/taxonomy.tsv", sep = "\t", header = TRUE)
  
  node_map = as.data.frame(full_taxonomy[,1])
  colnames(node_map) = c("Feature.ID")
  
  node_map = node_map %>%
    mutate(ASV_number = row_number())
  
  node_map_filtered <- node_map %>%
    filter(Feature.ID %in% rownames(subsampled_table))
  
  node_map_filtered_taxa = merge(node_map_filtered, sep.table_full, by.x = "Feature.ID", by.y = "row.names")
  
  ############### Make tidynet object
  
  tidy.net <- otu_net %>%
    as_tbl_graph()
  # Node size is proportional to the square root of the abundance of the corresponding OTU
  meanSize <- sqrt(rowMeans(subsampled_table)) %>%
    rescale(to = c(1,12))
  tidy.net %<>%
    activate(nodes) %>%
    mutate(nodeSize = unname(meanSize),
           asv_id = node_map_filtered_taxa$Feature.ID,
           name = node_map_filtered$ASV_number,
           counts = rowMeans(subsampled_table),
           taxon = node_map_filtered_taxa$Taxon)
  
  message("TidyNet object has been made.")
  
  ############### Delete isolated nodes
  
  rem_nodes <- which(degree(tidy.net)==0)
  tidy.net <- delete_vertices(tidy.net,rem_nodes)
  if(length(rem_nodes)!=0){
    final.tax <- final.tax[-c(rem_nodes)]
  }
  
  message("Isolated nodes have been deleted.")
  
  
  # Vertices
  vertices_df <- igraph::as_data_frame(tidy.net, what = "vertices")
  
  # Edges
  edges_df <- igraph::as_data_frame(tidy.net, what = "edges")
  
  # Graph-level attributes (as a named list)
  graph_attrs <- graph_attr(tidy.net)
  
  # Add number of edges to vertices table
  
  from_edge_count <- table(edges_df$from)
  to_edge_count <- table(edges_df$to)
  
  # Make everything 0 to get rid of NA values
  vertices_df$edge_count = 0
  
  
  for (node_name in vertices_df$name) {
    node_name <- trimws(node_name)  # clean whitespace
    
    from_count <- if (node_name %in% names(from_edge_count)) from_edge_count[[node_name]] else 0
    to_count   <- if (node_name %in% names(to_edge_count))   to_edge_count[[node_name]]   else 0
    
    vertices_df$edge_count[vertices_df$name == node_name] <- from_count + to_count
  }
  
  dir = file.path("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/", name)
  
  dir.create(dir, recursive = TRUE)
  
  write.table(vertices_df, file = file.path("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/", name, "vertices.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(edges_df, file = file.path("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/", name, "edges.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  graph_df <- data.frame(name = names(graph_attrs), value = unlist(graph_attrs))
  write.table(graph_df, file = file.path("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/", name, "graph_attributes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("TSVs of network attributes have been made.")
  
}

save.image("/lustre/isaac24/scratch/rsaho/projects/bsf_metaa/networks/meta_a_networks/animal_waste/animal_waste_network_gen.RData")
