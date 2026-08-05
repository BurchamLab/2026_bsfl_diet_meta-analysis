# Made by Reese Saho for the purpose of creating network figures based on the edges files created in the previous script.

invisible(gc()) # clear memory
rm(list = ls()) #clear environment

options(repos = c(CRAN = "https://cloud.r-project.org"))

packages = c("tidyverse", "tools", "igraph", "pals", "reshape2", "tidygraph", "scales", "ggraph", "ggdark", 
             "cowplot", 'dtplyr', 'googledrive', 'googlesheets4', 'haven', 'httr', 'ragg', 'rvest', 'xml2', 
             'Polychrome', 'GUniFrac', 'compositions', "reticulate", "ggplot2", "plotly", "RColorBrewer", "tidyr",
             "MESS", "parallel")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE) }
  
  library(pkg, character.only = TRUE)
}

message("Packages have been loaded.")

load("./sterile_plant_product/subsample/sterile_plant_product_network_gen.RData")

# Load Functions
bootstrap_edgedrop_modularity <- function(graph, n = 100, drop_fraction = 0.1) {
  mclapply(1:n, function(i) {
    # Sample edges to drop
    edges_to_drop <- sample(E(graph), size = floor(drop_fraction * ecount(graph)))
    subgraph <- delete_edges(graph, edges_to_drop)
    
    # Recalculate community and modularity
    mod_value <- modularity(cluster_fast_greedy(subgraph), directed = FALSE)
    return(mod_value)
  }, mc.cores = 4)
}

bootstrap_nodedrop_modularity <- function(graph, n = 100, drop_fraction = 0.1, method = "greedy") {
  # Arguments:
  # - graph: the igraph object to resample
  # - n: the number of resampling iterations (default: 10)
  # - drop_fraction: the fraction of nodes to remove each iteration (default: 10%)
  # - method: community detection method ("greedy", "louvain", "gn")
  
  # Ensure valid method
  if (!method %in% c("greedy", "louvain", "gn")) {
    stop("Method must be one of: 'greedy', 'louvain', or 'gn'")
  }
  
  # Function to select the community detection method
  detect_communities <- function(subgraph) {
    switch(method,
           greedy = cluster_fast_greedy(subgraph),
           louvain = cluster_louvain(subgraph),
           gn = cluster_edge_betweenness(subgraph, directed = FALSE))
  }
  
  # Store modularity values
  modularity_values <- mclapply(1:n, function(i) {
    # Randomly select vertices to remove
    vertices_to_drop <- sample(V(graph), size = floor(drop_fraction * vcount(graph)))
    subgraph <- delete_vertices(graph, vertices_to_drop)
    
    # Recalculate communities and modularity
    communities <- detect_communities(subgraph)
    modularity(communities, directed = FALSE)
  }, mc.cores = 4)
  
  # Return the list of modularity values
  return(modularity_values)
}

node_impact_analysis <- function(graph, method = "greedy") {
  # Arguments:
  # - graph: igraph object
  # - method: community detection method ("greedy", "louvain", "gn")
  
  # Original modularity and metrics
  original_modularity <- modularity(cluster_fast_greedy(graph), directed = FALSE)
  original_betweenness <- median(betweenness(graph, directed = FALSE))
  original_closeness <- median(closeness(graph))
  original_degree <- median(degree(graph))
  original_eigenvector <- median(eigen_centrality(graph)$vector, directed = FALSE)
  
  # Storage
  impact_df <- data.frame(
    Node = V(graph)$name,
    Modularity_Loss = NA,
    Betweenness_Loss = NA,
    Closeness_Loss = NA,
    Degree_Loss = NA,
    Eigenvector_Loss = NA
  )
  
  # Loop through each node
  for (v in V(graph)) {
    subgraph <- delete_vertices(graph, v)
    
    # Recalculate community detection and metrics
    new_modularity <- modularity(cluster_fast_greedy(subgraph), directed = FALSE)
    new_betweenness <- median(betweenness(subgraph, directed = FALSE))
    new_closeness <- median(closeness(subgraph))
    new_degree <- median(degree(subgraph))
    new_eigenvector <- median(eigen_centrality(subgraph)$vector, directed = FALSE)
    
    # Calculate absolute losses
    impact_df[impact_df$Node == V(graph)$name[v], "Modularity_Loss"] <- original_modularity - new_modularity
    impact_df[impact_df$Node == V(graph)$name[v], "Betweenness_Loss"] <- original_betweenness - new_betweenness
    impact_df[impact_df$Node == V(graph)$name[v], "Closeness_Loss"] <- original_closeness - new_closeness
    impact_df[impact_df$Node == V(graph)$name[v], "Degree_Loss"] <- original_degree - new_degree
    impact_df[impact_df$Node == V(graph)$name[v], "Eigenvector_Loss"] <- original_eigenvector - new_eigenvector
    
    # Calculate relative losses
    impact_df[impact_df$Node == V(graph)$name[v], "Modularity_Rel_Loss"] <- (original_modularity - new_modularity) / original_modularity
    impact_df[impact_df$Node == V(graph)$name[v], "Betweenness_Rel_Loss"] <- (original_betweenness - new_betweenness) / original_betweenness
    impact_df[impact_df$Node == V(graph)$name[v], "Closeness_Rel_Loss"] <- (original_closeness - new_closeness) / original_closeness
    impact_df[impact_df$Node == V(graph)$name[v], "Degree_Rel_Loss"] <- (original_degree - new_degree) / original_degree
    impact_df[impact_df$Node == V(graph)$name[v], "Eigenvector_Rel_Loss"] <- (original_eigenvector - new_eigenvector) / original_eigenvector
  }
  
  # Sort by impact
  impact_df <- impact_df[order(-impact_df$Modularity_Loss), ]
  rownames(impact_df) <- NULL
  return(impact_df)
}

resilience_test <- function(graph, impact_df, metric = NA, top_n = length(impact_df$Node), cumulative = TRUE) {
  
  # Ensure the metric is valid
  if (!(metric %in% colnames(impact_df) || metric == "random")) {
    stop("Invalid metric provided. Choose from: Modularity_Loss, Betweenness_Loss, Closeness_Loss, Degree_Loss, Eigenvector_Loss, Edge_Count_Loss, or random")
  }
  
  # Sort by the chosen metric and select the top nodes
  if (metric == "random") {
    # If random strategy, simply sample from the nodes
    top_nodes <- sample(impact_df$Node, top_n)
  } else {
    top_nodes <- impact_df %>%
      arrange(desc(.data[[metric]])) %>%
      head(top_n) %>%
      pull(Node)
  }
  
  resilience_metrics <- data.frame(
    Step = integer(),
    Components = integer(),
    Avg_Path_Length = numeric(),
    Clustering_Coeff = numeric()
  )
  
  if (cumulative) {
    # Cumulative removal: each step builds on the last
    for (i in 1:length(top_nodes)) {
      graph <- delete_vertices(graph, top_nodes[i])
      metrics <- data.frame(
        Step = i,
        Components = count_components(graph),
        Avg_Path_Length = mean_distance(graph, directed = FALSE),
        Clustering_Coeff = transitivity(graph, type = "global")
      )
      resilience_metrics <- rbind(resilience_metrics, metrics)
    }
  } else {
    # Independent removal: each node is removed from the original graph state
    for (i in 1:length(top_nodes)) {
      temp_graph <- delete.vertices(graph, top_nodes[i])
      metrics <- data.frame(
        Step = i,
        Components = count_components(temp_graph),
        Avg_Path_Length = mean_distance(temp_graph, directed = FALSE),
        Clustering_Coeff = transitivity(temp_graph, type = "global")
      )
      resilience_metrics <- rbind(resilience_metrics, metrics)
    }
  }
  return(resilience_metrics)
}

calculate_robustness_index <- function(graph, strategy = "degree", bootstrap = 100, method = "linear", remove_duplicates = TRUE) {
  # Calculate the initial number of edges and nodes
  initial_edges <- ecount(graph)
  initial_nodes <- vcount(graph)
  
  # If strategy is random, perform bootstrapping
  if (strategy == "random") {
    if (!is.numeric(bootstrap) || bootstrap <= 0) {
      stop("Invalid value for 'bootstrap'. It must be a positive number.")
    }
    robustness_indices <- c()
    edge_removal_lists <- list()
    giant_component_lists <- list()
    
    for (b in 1:bootstrap) {
      temp_graph <- graph
      edge_removal_fraction <- c()
      giant_component_fraction <- c()
      
      for (i in 1:initial_nodes) {
        if (vcount(temp_graph) > 0) {
          target_node <- sample(V(temp_graph), 1)
          temp_graph <- delete_vertices(temp_graph, target_node)
          
          edge_removal_fraction <- c(edge_removal_fraction, 1 - (ecount(temp_graph) / initial_edges))
          if (length(components(temp_graph)$csize) > 0) {
            giant_component_fraction <- c(giant_component_fraction, max(components(temp_graph)$csize) / initial_nodes)
          } else {
            giant_component_fraction <- c(giant_component_fraction, 0)
          }
        }
      }
      
      # Remove duplicates if specified
      if (remove_duplicates) {
        unique_indices <- !duplicated(edge_removal_fraction)
        edge_removal_fraction <- edge_removal_fraction[unique_indices]
        giant_component_fraction <- giant_component_fraction[unique_indices]
      }
      
      # Store the lists and calculate robustness for this bootstrap
      edge_removal_lists[[b]] <- edge_removal_fraction
      giant_component_lists[[b]] <- giant_component_fraction
      robustness_indices <- c(robustness_indices, auc(edge_removal_fraction, giant_component_fraction, type = method))
    }
    
    # Return the bootstrapped results
    return(list(
      Strategy = strategy,
      Robustness_Index = robustness_indices,
      Edge_Removal_Fraction = edge_removal_lists,
      Giant_Component_Fraction = giant_component_lists
    ))
  }
  
  # If not random, proceed with original method
  edge_removal_fraction <- c()
  giant_component_fraction <- c()
  
  for (i in 1:initial_nodes) {
    if (vcount(graph) > 0) {
      if (strategy == "degree") {
        target_node <- which.max(degree(graph))
      } else if (strategy == "betweenness") {
        target_node <- which.max(betweenness(graph, directed = FALSE))
      } else if (strategy == "closeness") {
        target_node <- which.max(closeness(graph))
      } else if (strategy == "eigenvector") {
        target_node <- which.max(eigen_centrality(graph)$vector, directed = FALSE)
      } else {
        stop("Invalid strategy. Choose from: degree, betweenness, closeness, eigenvector, random")
      }
      
      graph <- delete_vertices(graph, target_node)
      
      # Record the fractions
      edge_removal_fraction <- c(edge_removal_fraction, 1 - (ecount(graph) / initial_edges))
      if (length(components(graph)$csize) > 0) {
        giant_component_fraction <- c(giant_component_fraction, max(components(graph)$csize) / initial_nodes)
      } else {
        giant_component_fraction <- c(giant_component_fraction, 0)
      }
    }
  }
  
  # Remove duplicates if specified
  if (remove_duplicates) {
    unique_indices <- !duplicated(edge_removal_fraction)
    edge_removal_fraction <- edge_removal_fraction[unique_indices]
    giant_component_fraction <- giant_component_fraction[unique_indices]
  }
  
  # Calculate the AUC with the trapezoidal rule (linear)
  robustness_index <- auc(edge_removal_fraction, giant_component_fraction, type = method)
  
  # Return the index and the breakdown
  return(list(
    Strategy = strategy,
    Robustness_Index = robustness_index,
    Edge_Removal_Fraction = edge_removal_fraction,
    Giant_Component_Fraction = giant_component_fraction
  ))
}

interactive_impact_plots <- function(impact_df) {
  # Combined Plots
  impact_long <- impact_df %>%
    tidyr::pivot_longer(cols = c(Modularity_Loss, Betweenness_Loss, Closeness_Loss, Degree_Loss, Eigenvector_Loss),
                        names_to = "Metric", values_to = "Loss")
  impact_rel_long <- impact_df %>%
    tidyr::pivot_longer(cols = c(Modularity_Rel_Loss, Betweenness_Rel_Loss, Closeness_Rel_Loss, Degree_Rel_Loss, Eigenvector_Rel_Loss),
                        names_to = "Metric", values_to = "Relative_Loss")
  
  p1 <- plot_ly(data = impact_long, x = ~Node, y = ~Loss, color = ~Metric, type = 'bar') %>%
    layout(title = 'Impact of Node Removal Across Metrics',
           xaxis = list(title = 'Node'),
           yaxis = list(title = 'Loss'))
  
  p2 <- plot_ly(data = impact_rel_long, x = ~Node, y = ~Relative_Loss, color = ~Metric, type = 'bar') %>%
    layout(title = 'Relative Impact of Node Removal Across Metrics',
           xaxis = list(title = 'Node'),
           yaxis = list(title = 'Relative Loss'))
  
  # Return the list of plots
  return(list(Combined = p1, Combined_Relative = p2))
}


# Set the number of iterations that was conducted

iterations = 100

# Make an edge list file that merges all of the made edge files

base_path = "./sterile_plant_product/subsample"

edge_files = list.files(path = base_path, pattern = ".*edges\\.tsv$", 
                       recursive = TRUE, 
                       full.names = TRUE)

edge_data_list <- lapply(edge_files, read.delim)

message("List of edge files has been made.")

edge_file_full = bind_rows(edge_data_list)

# Make a new column that shows interactions as "ASV 1" - "ASV 2"

edge_file_full$interactions = paste(edge_file_full$from, edge_file_full$to, sep = "-")

# Make a new column that calculates how frequently an interaction was present in an iteration

edge_file_full = edge_file_full %>%
  group_by(from, to) %>%
  summarise(interaction_count = n() / iterations, .groups = "drop")

edge_file_filtered = edge_file_full %>%
  filter(interaction_count >= 0.5)


sterile_plant_product_graph = graph_from_data_frame(edge_file_filtered, directed = FALSE)

node_map_filtered <- node_map %>%
  filter(ASV_number %in% V(sterile_plant_product_graph)$name)

V(sterile_plant_product_graph)$asv_id = node_map_filtered$Feature.ID

# 1. Community Detection (Greedy)
comm_greedy <- cluster_fast_greedy(sterile_plant_product_graph)

# 2. Modularity Calculation
mod_greedy <- modularity(comm_greedy, directed = FALSE)

cat("Greedy Modularity - Sterile Plant Product Diet:", mod_greedy)

# 3. Centrality Analysis
# Calculate betweenness centrality to identify hubs
betweenness_values <- betweenness(sterile_plant_product_graph, normalized = TRUE, directed = FALSE)
V(sterile_plant_product_graph)$betweenness <- betweenness_values

# Extract top 5 hub nodes
top_hubs <- V(sterile_plant_product_graph)[order(-V(sterile_plant_product_graph)$betweenness)][1:5]$name

# Calculate betweenness centrality to identify edges
edge_betweenness_values <- edge_betweenness(sterile_plant_product_graph, directed = FALSE)
E(sterile_plant_product_graph)$edge_betweenness <- edge_betweenness_values

# Extract top 5 central edges
top_edges <- E(sterile_plant_product_graph)[order(-E(sterile_plant_product_graph)$edge_betweenness)][1:5]
edge_names <- sapply(top_edges, function(edge) paste(ends(sterile_plant_product_graph, edge), collapse = "--"))

cat("Central Nodes in Sterile Plant Product Diet:", top_hubs, "\n")
cat("Central Edges in Sterile Plant Product Diet:", edge_names, "\n")

# 4. Node Impact Analysis
impact <- node_impact_analysis(sterile_plant_product_graph)
impact_plots <- interactive_impact_plots(impact)
impact_plots$Combined
impact_plots$Combined_Relative

# 5. Resilience Test based on Impactful Nodes
cumulative_resilience <- resilience_test(sterile_plant_product_graph, impact, metric = "Modularity_Loss", cumulative = TRUE)
cumulative_resilience_random <- resilience_test(sterile_plant_product_graph, impact, metric = "random", cumulative = TRUE)


# Plot the number of components over time
ggplot(cumulative_resilience, aes(x = Step, y = Components)) +
  geom_line(color = "blue") +
  labs(title = "Number of Components During Node Removal", x = "Step (Node Removal)", y = "Number of Components") +
  theme_minimal()
ggsave(filename = "./sterile_plant_product/sterile_plant_product_components.pdf", device = "pdf", height = 8, width = 10, units = "in")

ggplot(cumulative_resilience_random, aes(x = Step, y = Components)) +
  geom_line(color = "blue") +
  labs(title = "Number of Components During Node Removal", x = "Step (Node Removal)", y = "Number of Components") +
  theme_minimal()
ggsave(filename = "./sterile_plant_product/sterile_plant_product_components_random.pdf", device = "pdf", height = 8, width = 10, units = "in")


# Plot the average path length over time
ggplot(cumulative_resilience, aes(x = Step, y = Avg_Path_Length)) +
  geom_line(color = "red") +
  labs(title = "Average Path Length During Node Removal", x = "Step (Node Removal)", y = "Average Path Length") +
  theme_minimal()
ggsave(filename = "./sterile_plant_product/sterile_plant_product_avg_path_length.pdf", device = "pdf", height = 8, width = 10, units = "in")


ggplot(cumulative_resilience_random, aes(x = Step, y = Avg_Path_Length)) +
  geom_line(color = "red") +
  labs(title = "Average Path Length During Node Removal", x = "Step (Node Removal)", y = "Average Path Length") +
  theme_minimal()
ggsave(filename = "./sterile_plant_product/sterile_plant_product_avg_path_length_random.pdf", device = "pdf", height = 8, width = 10, units = "in")


# Plot the clustering coefficient over time
ggplot(cumulative_resilience, aes(x = Step, y = Clustering_Coeff)) +
  geom_line(color = "green") +
  labs(title = "Global Clustering Coefficient During Node Removal", x = "Step (Node Removal)", y = "Clustering Coefficient") +
  theme_minimal()
ggsave(filename = "./sterile_plant_product/sterile_plant_product_global_clust_coef.pdf", device = "pdf", height = 8, width = 10, units = "in")

ggplot(cumulative_resilience_random, aes(x = Step, y = Clustering_Coeff)) +
  geom_line(color = "green") +
  labs(title = "Global Clustering Coefficient During Node Removal", x = "Step (Node Removal)", y = "Clustering Coefficient") +
  theme_minimal()
ggsave(filename = "./sterile_plant_product/sterile_plant_product_global_clust_coef_random.pdf", device = "pdf", height = 8, width = 10, units = "in")


# 6. Robustness Test
# Run the robustness calculation using strategy to remove nodes in oder of highest metric (e.g., 'degree' removals node with highest edge count then node with 2nd highest etc.)
robustness_results_degree <- calculate_robustness_index(sterile_plant_product_graph, strategy = "degree", method = "linear")
robustness_results_betweenness <- calculate_robustness_index(sterile_plant_product_graph, strategy = "betweenness", method = "linear")

# Run the robustness calculation using random strategy so that a node is chosen at random to remove, this is bootstrapped to tests random variability
robustness_rand_results <- calculate_robustness_index(sterile_plant_product_graph, strategy = "random", bootstrap = 100, method = "linear")

# View the robustness index
cat("Robustness AUC - Sterile Plant Product Diet (degree):", robustness_results_degree$Robustness_Index, "\n")
cat("Robustness AUC- Sterile Plant Product Diet (betweenness):", robustness_results_betweenness$Robustness_Index, "\n")

# Convert to DataFrame for plotting
plot_data <- data.frame(
  Step = 1:length(robustness_results_degree$Giant_Component_Fraction),
  Giant_Component = robustness_results_degree$Giant_Component_Fraction,
  Edge_Removal = robustness_results_degree$Edge_Removal_Fraction
)

# Plot with shaded area and AUC in the title
ggplot(plot_data, aes(x = Edge_Removal, y = Giant_Component)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_area(fill = "lightblue", alpha = 0.3) +
  labs(
    title = paste0("Robustness Analysis (AUC = ", round(robustness_results_degree$Robustness_Index, 4), ")"),
    x = "Edge Removal Fraction",
    y = "Giant Component Fraction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
ggsave(filename = "./sterile_plant_product/sterile_plant_product_robustness_analysis.pdf", device = "pdf", height = 8, width = 10, units = "in")

# Plot Random
# Combine all bootstrap results into a data frame
random_trials <- robustness_rand_results$Giant_Component_Fraction
edge_removal_trials <- robustness_rand_results$Edge_Removal_Fraction

# Create a long-format data frame for plotting
bootstrap_df <- do.call(rbind, lapply(seq_along(random_trials), function(i) {
  data.frame(
    Trial = paste("Bootstrap", i),
    Step = 1:length(random_trials[[i]]),
    Giant_Component = random_trials[[i]],
    Edge_Removal = edge_removal_trials[[i]]
  )
}))


# Calculate the mean decay across all trials
mean_df <- bootstrap_df %>%
  group_by(Step) %>%
  summarize(
    Mean_Giant_Component = mean(Giant_Component),
    Mean_Edge_Removal = mean(Edge_Removal)
  )

# Calculate Mean AUC
mean_auc <- mean(robustness_rand_results$Robustness_Index, na.rm = TRUE)
cat("Mean Random Robustness AUC:", mean_auc, "\n")

# Plot all bootstrap lines and the average with shaded area
ggplot() +
  # Plot all bootstrapped lines in light gray
  geom_line(data = bootstrap_df, aes(x = Edge_Removal, y = Giant_Component, group = Trial), 
            color = "gray", alpha = 0.6) +
  # Plot the filled area under the mean curve
  geom_ribbon(data = mean_df, aes(x = Mean_Edge_Removal, ymin = 0, ymax = Mean_Giant_Component), 
              fill = "lightblue", alpha = 0.3) +
  # Plot the mean line on top
  geom_line(data = mean_df, aes(x = Mean_Edge_Removal, y = Mean_Giant_Component), 
            color = "blue", size = 1.2, linetype = "solid") +
  
  labs(
    title = paste0("Random Removal Robustness Analysis Sterile Plant Product Diet (Mean AUC: ", round(mean_auc, 4), ")"),
    subtitle = "Individual Trials (light gray), Mean Curve (blue), Shaded AUC (light blue)",
    x = "Edge Removal Fraction",
    y = "Giant Component Fraction"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
ggsave(filename = "./sterile_plant_product/sterile_plant_product_waste_random_removal_robustness.pdf", device = "pdf", width = 10, height = 8, units = "in")

# 6. Visualization
par(mfrow = c(2, 2))

# Generate community colors
community_colors <- rainbow(length(unique(membership(comm_greedy))))

# Get the membership of each node
memberships <- membership(comm_greedy)

# Map each node to its corresponding community color
V(sterile_plant_product_graph)$color <- community_colors[memberships]

# Combine and collapse all of the generated feature tables

base_path = "./sterile_plant_product/"
# base_path = "/Users/burchamlab/Library/CloudStorage/OneDrive-UniversityofTennessee/UTK/Research/Code/bsf_metaa/SparCC-master/sterile_plant_product/subsample/"

feature_table_files = list.files(path = base_path, pattern = ".*_abun_[0-9]+\\.txt$", 
                       recursive = TRUE, 
                       full.names = TRUE)

feature_table_list <- lapply(feature_table_files, read.delim, row.names = 1)

message("List of correlation files has been made.")

names(feature_table_list) <- stringr::str_extract(feature_table_files, "(?<=abun_)[0-9]+")
  

# Step 1: Convert rownames to column "ASV"
feature_table_list_named <- lapply(feature_table_list, function(df) {
  df %>% rownames_to_column("ASV")
})

# Step 2: Bind all rows
feature_table_merged <- reduce(feature_table_list_named, .init = tibble(ASV = character()), full_join, by = "ASV")

# Step 3: Remove duplicate ASVs (keep first occurrence)
feature_table_merged <- column_to_rownames(feature_table_merged, "ASV")

# Transpose the merged table so samples are rows
transposed <- t(feature_table_merged)

# Remove duplicate rows (i.e., samples with identical ASV abundances)
unique_transposed <- transposed[!duplicated(as.data.frame(transposed)), ]

# Transpose back so ASVs are rows again
feature_table_deduped <- t(unique_transposed)

# Remove common suffixes like .x, .y, .1, .2, etc.
clean_names <- gsub("\\.[xy]|\\.\\d+$", "", colnames(feature_table_deduped))

# Apply the cleaned names back
colnames(feature_table_deduped) <- clean_names


# Make a graph that colors by study

t_fft = as.data.frame(t(feature_table_deduped))

t_fft$SampleID <- rownames(t_fft)
metadata$SampleID <- rownames(metadata)

merged_table_study <- merge(metadata, t_fft, by = "SampleID")

merged_table_study[7:ncol(merged_table_study)] <- lapply(
  merged_table_study[7:ncol(merged_table_study)],
  as.numeric
)

collasped_ft_study = merged_table_study %>%
  group_by(study_id) %>%
  summarise(across(6:(ncol(merged_table_study)-1), \(x) mean(x, na.rm = TRUE)))

t_cft_study = as.data.frame(t(collasped_ft_study), header = TRUE)

colnames(t_cft_study) <- as.character(t_cft_study[1, ])
t_cft_study <- t_cft_study[-1, , drop = FALSE]

matched_t_cft_study <- t_cft_study[match(V(sterile_plant_product_graph)$asv_id, rownames(t_cft_study)), , drop = FALSE]
V(sterile_plant_product_graph)$study_pie <- lapply(asplit(matched_t_cft_study, 1), as.numeric)

study_colors = createPalette(ncol(matched_t_cft_study), seedcolors = (c("#FF0000", "#0000FF", "#00FF00")))
names(study_colors) = colnames(t_cft_study)
study_colors <- study_colors[!is.na(names(study_colors))]  # Drop unnamed colors


set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_study_network.pdf", width = 10, height = 8)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.shape = "pie",
     vertex.pie = V(sterile_plant_product_graph)$study_pie,
     vertex.pie.color = list(study_colors),
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - Study Level"
)

legend("topright", legend = names(study_colors), fill = study_colors, border = "black", cex = 0.8)

dev.off()

############### Make a diet-based graph with pie chart nodes

merged_table_diet <- merge(metadata, t_fft, by = "SampleID")

merged_table_diet[8:ncol(merged_table_diet)] <- lapply(
  merged_table_diet[8:ncol(merged_table_diet)],
  as.numeric
)

collasped_ft_diet = merged_table_diet %>%
  group_by(sterile_plant_product_type) %>%
  summarise(across(8:(ncol(merged_table_diet) - 1), \(x) mean(x, na.rm = TRUE)))

t_cft_diet = as.data.frame(t(collasped_ft_diet), header = TRUE)

colnames(t_cft_diet) <- as.character(t_cft_diet[1, ])
t_cft_diet <- t_cft_diet[-1, , drop = FALSE]

matched_t_cft_diet <- t_cft_diet[match(V(sterile_plant_product_graph)$asv_id, rownames(t_cft_diet)), , drop = FALSE]
V(sterile_plant_product_graph)$diet_pie <- lapply(asplit(matched_t_cft_diet, 1), as.numeric)

diet_colors = createPalette(ncol(matched_t_cft_diet), seedcolors = (c("#8B0000", "#00008B", "#006400")))
names(diet_colors) = colnames(t_cft_diet)
diet_colors <- diet_colors[!is.na(names(diet_colors))]  # Drop unnamed colors


set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_diet_network.pdf", width = 10, height = 8)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.shape = "pie",
     vertex.pie = V(sterile_plant_product_graph)$diet_pie,
     vertex.pie.color = list(diet_colors),
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - Diet Level"
)

legend("topright", legend = names(diet_colors), fill = diet_colors, border = "black", cex = 0.8)

dev.off()


# Plot with matching node and community colors

set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_cluster_network.pdf", width = 8, height = 8)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(comm_greedy, sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.size = ifelse(V(sterile_plant_product_graph)$name %in% top_hubs, 10, 3),  # Hubs larger
     vertex.label = ifelse(V(sterile_plant_product_graph)$name %in% top_hubs, V(sterile_plant_product_graph)$name, NA),  # Only label hubs
     vertex.label.cex = 0.8,  # Change label size
     vertex.label.dist = 0,  # Adjust label distance from the node
     vertex.label.color = "black",  # Label color for visibility
     col = V(sterile_plant_product_graph)$color,  # Explicitly set node colors
     mark.groups = communities(comm_greedy),  # Mark the communities
     mark.col = adjustcolor(community_colors, alpha.f = 0.4),  # Community boundary color with transparency
     mark.border = community_colors,  # Border of community circles
     main = "sterile_plant_product Waste Network - Clusters & Top Hubs Highlighted",
     edge.color = "gray")

dev.off()

############### Make a graph that shows labels by ASVs

V(sterile_plant_product_graph)$asv_name <- V(tidy.net)$asv_name[match(V(sterile_plant_product_graph)$name, V(tidy.net)$name)]

ASV_names <- unique(V(sterile_plant_product_graph)$asv_name)

# Separate "Unknown" and other phyla
known_ASVs <- sort(ASV_names[ASV_names != "Unknown"])
ASV_list <- c(known_ASVs, "Unknown")  # put "Unknown" last


ASV_colors = createPalette(length(known_ASVs), seedcolors = (c("#8B0000", "#00008B", "#006400")))
names(ASV_colors) = known_ASVs

ASV_colors["Unknown"] <- "#555555"

V(sterile_plant_product_graph)$ASV_color <- ASV_colors[V(sterile_plant_product_graph)$asv_name]

legend_labels <- ASV_list

set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_ASV_network.pdf", width = 18, height = 10)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.color = V(sterile_plant_product_graph)$ASV_color,
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - ASV Level"
)

legend("topright", legend = legend_labels, fill = ASV_colors, border = "black", cex = 0.45)

dev.off()

############### Make a graph that colors by genus


V(sterile_plant_product_graph)$genus_name = tax$Genus[match(V(sterile_plant_product_graph)$asv_id, rownames(tax))]

genus_names <- unique(V(sterile_plant_product_graph)$genus_name)

# Separate "Unknown" and other phyla
known_genera <- sort(genus_names[genus_names != "Unknown"])
genus_list <- c(known_genera, "Unknown")  # put "Unknown" last


Genus_colors = createPalette(length(known_genera), seedcolors = (c("#8B0000", "#00008B", "#006400")))
names(Genus_colors) = known_genera

Genus_colors["Unknown"] <- "#555555"

V(sterile_plant_product_graph)$genus_color <- Genus_colors[V(sterile_plant_product_graph)$genus_name]

legend_labels <- genus_list

set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_genus_network.pdf", width = 14, height = 10)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.color = V(sterile_plant_product_graph)$genus_color,
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - Genus Level"
)

legend("topright", legend = legend_labels, fill = Genus_colors, border = "black", cex = 0.7)

dev.off()

############### Make a graph that colors by family


V(sterile_plant_product_graph)$family_name = tax$Family[match(V(sterile_plant_product_graph)$asv_id, rownames(tax))]

family_names <- unique(V(sterile_plant_product_graph)$family_name)

# Separate "Unknown" and other phyla
known_families <- sort(family_names[family_names != "Unknown"])
family_list <- c(known_families, "Unknown")  # put "Unknown" last


Family_colors = createPalette(length(known_families), seedcolors = (c("#8B0000", "#00008B", "#006400")))
names(Family_colors) = known_families

Family_colors["Unknown"] <- "#555555"

V(sterile_plant_product_graph)$family_color <- Family_colors[V(sterile_plant_product_graph)$family_name]

legend_labels <- family_list

set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_family_network.pdf", width = 14, height = 10)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.color = V(sterile_plant_product_graph)$family_color,
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - Family Level"
)

legend("topright", legend = legend_labels, fill = Family_colors, border = "black", cex = 0.7)

dev.off()

############### Make a graph that colors by order


V(sterile_plant_product_graph)$order_name = tax$Order[match(V(sterile_plant_product_graph)$asv_id, rownames(tax))]

order_names <- unique(V(sterile_plant_product_graph)$order_name)

# Separate "Unknown" and other phyla
known_orders <- sort(order_names[order_names != "Unknown"])
order_list <- c(known_orders, "Unknown")  # put "Unknown" last


Order_colors = createPalette(length(known_orders), seedcolors = (c("#8B0000", "#00008B", "#006400")))
names(Order_colors) = known_orders

Order_colors["Unknown"] <- "#555555"

V(sterile_plant_product_graph)$order_color <- Order_colors[V(sterile_plant_product_graph)$order_name]

legend_labels <- order_list

set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_order_network.pdf", width = 14, height = 10)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.color = V(sterile_plant_product_graph)$order_color,
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - Order Level"
)

legend("topright", legend = legend_labels, fill = Order_colors, border = "black", cex = 0.7)

dev.off()

############### Make a graph that colors by class


V(sterile_plant_product_graph)$class_name = tax$Class[match(V(sterile_plant_product_graph)$asv_id, rownames(tax))]

class_names <- unique(V(sterile_plant_product_graph)$class_name)

# Separate "Unknown" and other phyla
known_classes <- sort(class_names[class_names != "Unknown"])
class_list <- c(known_classes, "Unknown")  # put "Unknown" last


Class_colors = createPalette(length(known_classes), seedcolors = (c("#8B0000", "#00008B", "#006400")))
names(Class_colors) = known_classes

Class_colors["Unknown"] <- "#555555"

V(sterile_plant_product_graph)$class_color <- Class_colors[V(sterile_plant_product_graph)$class_name]

legend_labels <- class_list

set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_class_network.pdf", width = 14, height = 10)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.color = V(sterile_plant_product_graph)$class_color,
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - Class Level"
)

legend("topright", legend = legend_labels, fill = Family_colors, border = "black", cex = 0.7)

dev.off()

############### Make a graph that colors by Phylum


V(sterile_plant_product_graph)$phylum_name = tax$Phylum[match(V(sterile_plant_product_graph)$asv_id, rownames(tax))]

phylum_names <- unique(V(sterile_plant_product_graph)$phylum_name)

# Separate "Unknown" and other phyla
known_phyla <- sort(phylum_names[phylum_names != "Unknown"])
phylum_list <- c(known_phyla, "Unknown")  # put "Unknown" last


Phylum_colors = createPalette(length(known_phyla), seedcolors = (c("#8B0000", "#00008B", "#006400")))
names(Phylum_colors) = known_phyla

Phylum_colors["Unknown"] <- "#555555"

V(sterile_plant_product_graph)$phylum_color <- Phylum_colors[V(sterile_plant_product_graph)$phylum_name]

legend_labels <- phylum_list

set.seed(76)
pdf("./sterile_plant_product/sterile_plant_product_phylum_network.pdf", width = 14, height = 10)

layout = layout_nicely(sterile_plant_product_graph)
layout_expanded = layout * 5

plot(sterile_plant_product_graph,
     layout = layout_expanded,
     vertex.color = V(sterile_plant_product_graph)$phylum_color,
     vertex.size = 4,
     vertex.label.color = "white",
     edge.color = "darkgrey",
     vertex.label.cex = 0.35,
     main = "Sterile Plant Product Network - Phylum Level"
)

legend("topright", legend = legend_labels, fill = Phylum_colors, border = "black", cex = 0.7)

dev.off()

# Make a master edges and vertices file from the overall network

# Vertices
vertices_df <- igraph::as_data_frame(sterile_plant_product_graph, what = "vertices")
non_list_cols <- sapply(vertices_df, function(col) !is.list(col)) # Removes study and diet pie charts, can't save them because they are a list
vertices_df_clean <- vertices_df[, non_list_cols]
write.table(vertices_df_clean, file = file.path("./sterile_plant_product/sterile_plant_product_vertices.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# Edges
edges_df <- igraph::as_data_frame(sterile_plant_product_graph, what = "edges")
write.table(edges_df, file = file.path("./sterile_plant_product/sterile_plant_product_edges.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

