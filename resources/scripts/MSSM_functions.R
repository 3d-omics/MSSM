# print_median_iqr(
print_median_iqr <- function(x, name) {
  med <- round(median(x), 2)
  iqr_vals <- round(quantile(x, probs = c(0.25, 0.75)), 2)
  cat(name, "median =", med, "IQR =", iqr_vals[1], "-", iqr_vals[2], "\n")
}


# calculate_alpha_diversity()
calculate_alpha_diversity <- function(input_data, metadata_name, tree_name) {
  input_data_matrix <- input_data %>% column_to_rownames(var = "genome")
  colsum <- input_data_matrix %>%
    summarise(across(where(is.numeric), sum)) %>%
    pivot_longer(cols = everything(), names_to = "microsample", values_to = "cov_filtering") %>%
    mutate(cov_filtering = ifelse(cov_filtering > 0, "Retained By Filtering", "Excluded By Filtering"))
  # Define diversity metrics and calculate
  diversity_metrics <- list(
    richness = hilldiv(input_data_matrix, q = 0),
    neutral = hilldiv(input_data_matrix, q = 1),
    phylogenetic = hilldiv(input_data_matrix, q = 1, tree = tree_name)
  )
  # Process metrics into a single data frame
  alpha_diversity <- lapply(names(diversity_metrics), function(metric) {
    diversity_metrics[[metric]] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = "microsample") %>%
      rename(!!sym(metric) := 2) # rename metric column
  }) %>%
    reduce(full_join, by = "microsample") %>%
    right_join(metadata_name, by = "microsample") %>% # Merge with final stats
    left_join(colsum, by = "microsample")
  return(alpha_diversity)
}

# MAG_tree_plot()
### Generate the phylum color heatmap
MAG_tree_plot <- function(genome_metadata, tree) {
  phylum_heatmap <- genome_metadata %>%
    arrange(match(genome, tree$tip.label)) %>%
    select(genome,phylum) %>%
    mutate(phylum = factor(phylum, levels = unique(phylum))) %>%
    column_to_rownames(var = "genome")
  
  ### Generate the order color heatmap
  order_heatmap <- genome_metadata %>%
    arrange(match(genome, tree$tip.label)) %>%
    select(genome, order) %>%
    mutate(order = factor(order, levels = unique(order))) %>%
    column_to_rownames(var = "genome")
  
  ### Generate the basal tree
  circular_tree <- force.ultrametric(tree, method="extend") %>% # extend to ultrametric for the sake of visualisation
    ggtree(., layout="fan", open.angle=10, size=0.5)
  
  ### Add the phylum ring
  circular_tree <- gheatmap(circular_tree, 
                            phylum_heatmap, 
                            offset = 1.2, width = 0.2, colnames = FALSE) +
    scale_fill_manual(values = phylum_colors, name = "Phylum") +
    geom_tiplab2(size = 1, hjust = -0.1) +
    theme(plot.margin = margin(0, 0, 0, 0), 
          panel.margin = margin(0, 0, 0, 0))
    circular_tree <- circular_tree + new_scale_fill()
  
  ### Add order ring
  circular_tree <- gheatmap(circular_tree, 
                            order_heatmap, 
                            offset = 0.55, width = 0.3, colnames = FALSE) +
    scale_fill_manual(values = order_colors, name = "Order") +
    theme(legend.position = "bottom", legend.box = "vertical") 
    circular_tree <- circular_tree + new_scale_fill()
  
  ### Add completeness ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_gradient(low = "#d9727f", high = "#465c7a", 
                        name = "Completeness") +
    geom_fruit(data=genome_metadata, geom=geom_bar, 
               mapping = aes(x=completeness, y=genome, fill=completeness),
               offset = 1.00, orientation="y", stat="identity")
  
  ### Add contamination ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_gradient(low = "#465c7a", high = "#d9727f",
                        name = "Contamination %") +
    geom_fruit(data=genome_metadata, geom=geom_bar,
               mapping = aes(x=contamination, y=genome, fill=contamination),
               offset = 0.05, orientation="y", stat="identity")
  
  ### Add genome size ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_gradient(low = "#465c7a", high = "#fcbb6d", 
                        name = "Genome size") +
    geom_fruit(data = genome_metadata, geom = geom_bar,
               mapping = aes(x = length, y = genome, fill = length),
               offset = 0.01, orientation = "y", stat = "identity")
  
  ### Add circularised genome ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_manual(values = c("1" = "#fcbb6d", "0" = "#465c7a"), 
                      name = 'Circularized genomes') +
    geom_fruit(data = genome_metadata, geom = geom_bar,
               mapping = aes(x = 1, y = genome, fill = circularity),
               offset = 0.02, orientation = "y", stat = "identity")
  
  # Add text
  circular_tree <-  circular_tree +
    annotate('text', x=3.1, y=0, label='           Order', size=3.5) +
    annotate('text', x=3.55, y=0, label='             Phylum', size=3.5) +
    annotate('text', x=4.1, y=0, label='                         Completeness', size=3.5) +
    annotate('text', x=4.5, y=0, label='                         Contamination', size=3.5) +
    annotate('text', x=5.0, y=0, label='                       Genome size', size=3.5) +
    annotate('text', x=5.5, y=0, label='                      Circularised', size=3.5)
  
  # Plot circular tree
  p <-circular_tree %>% open_tree(30) %>% rotate_tree(90) 
  return(p)
}


# remove_samples_or_taxa()
## Input: A wide df with samples as rows and taxa as columns
## Result: Remove taxa that are in too few samples, or samples that have too few taxa
remove_samples_or_taxa <- function(df, min_samples_per_taxon, min_taxa_per_sample){
  # Store & print original df dimensions
  original_rows <- nrow(df)
  original_cols <- ncol(df)
  cat("Initial df: Rows (samples):", original_rows, ", Columns (taxa):", original_cols, "\n")
  
  # Remove taxa that are in less than x samples
  df <- df %>% select(where(~ sum(. != 0) > min_samples_per_taxon))
  
  # Remove samples that contain less than x taxa
  df <- df %>% filter(rowSums(. != 0) > min_taxa_per_sample)
  
  df < df %>% select(where(~ any(. != 0)) | where(is.character) | where(is.factor))  
  
  removed_rows <- original_rows - nrow(df)
  removed_cols <- original_cols - ncol(df)
  cat("Removed: Rows (samples):", removed_rows, ", Columns (taxa):", removed_cols, "\n")
  
  cat("Resulting df: Rows (samples):", nrow(df), ", Columns (taxa):", ncol(df), "\n")
  
  return(df)
}


# clr_transform
clr_transform <- function(x) {
  log(x) - mean(log(x), na.rm = TRUE)
}

# perform_pca()
perform_pca <- function(df, zero_method = "GBM", z_delete = TRUE) {
  # Store original dimensions
  original_rows <- nrow(df)
  original_cols <- ncol(df)
  
  # 1. Zero replacement
  if (any(df == 0)) { # I think cmultRepl already does that
    print("Zeros found")
    df <- cmultRepl(df,
                    method = zero_method, output = "prop",
                    z.warning = 0.8, z.delete = z_delete
    )
    df <- df * 100
  }
  
  # Print removed rows and columns
  removed_rows <- original_rows - nrow(df)
  removed_cols <- original_cols - ncol(df)
  cat("Rows (samples) removed after zero replacement:", removed_rows, "\n")
  cat("Columns (taxa) removed after zero replacement:", removed_cols, "\n")
  
  # Geometric mean function
  geometric_mean <- function(x) {
    # Use log to avoid underflow
    exp(mean(log(x), na.rm = TRUE))
  }
  
  # 2. Calculate geometric mean of the parts (taxa) of the data set.
  taxa_geometric_means <- apply(df, 2, geometric_mean)
  
  # 3. Center data
  df_centered <- sweep(df, 2, taxa_geometric_means, FUN = "/")
  
  df_centered <- as.matrix(df_centered)
  
  # Compute the Variation Matrix
  variation_matrix <- outer(
    1:ncol(df_centered), 1:ncol(df_centered),
    Vectorize(function(i, j) var(log(df_centered[, i] / df_centered[, j]), na.rm = TRUE))
  )
  
  # Calculate Total Variance
  D <- ncol(df_centered) # Number of taxa (columns)
  totvar <- (1 / (2 * D)) * sum(variation_matrix, na.rm = TRUE)
  
  # 4. Scale data
  power_exponent <- 1 / sqrt(totvar)
  df_scaled <- df_centered^power_exponent
  
  # CLR transform data
  df_clr <- as.data.frame(t(apply(df_scaled, 1, clr_transform)))
  
  df_clr_dist <- as.data.frame(t(apply(df, 1, clr_transform)))
  
  # Perform PCA on zero replaced, centered, scaled, and CLR transformed df
  pca_result <- prcomp(df_clr, center = FALSE, scale. = FALSE)
  
  return(list(
    df_clr = df_clr,
    df_clr_dist = df_clr_dist,
    pca_result = pca_result
  ))
}

# plot_pca() 
### Use the 'pca_result' df produced from 'perform_pca' function to make the PCA plot
plot_pca <- function(df, 
                     samples_color_metadata, samples_shape_metadata, 
                     samples_color_value, loadings_color_metadata, 
                     loadings_color_value, loadings_taxon_level,
                     sample_metadata, genome_metadata, order_colors, 
                     custom_ggplot_theme, scaling_factor_value = 1.5, 
                     loadings_number = 10) {
  
  # Extract scores from PCA results
  scores <- rownames_to_column(as.data.frame(df$x), var = "microsample")
  scores <- left_join(scores, sample_metadata, by = join_by(microsample == microsample))
  
  # Calculate limits for x and y axes
  x_limit <- max(abs(scores$PC1))
  y_limit <- max(abs(scores$PC2))
  
  # Calculate variance explained by each PC (principal component) & create labels for plot
  variance_explained <- (df$sdev^2) / sum(df$sdev^2) * 100
  pc1_label <- paste0("PC1: ", round(variance_explained[1], 2), "% variance explained")
  pc2_label <- paste0("PC2: ", round(variance_explained[2], 2), "% variance explained")
  
  # Set a scaling factor for loadings (arrows)
  scaling_factor <- scaling_factor_value
  
  # Extract and scale the loadings
  loadings <- df$rotation[, 1:2] %>%
    as.data.frame() %>%
    mutate(genome = rownames(.)) %>%  
    mutate(PC1 = PC1 * scaling_factor,
           PC2 = PC2 * scaling_factor) %>%
    left_join(genome_metadata, by = join_by(genome == genome)) %>%  
    mutate(abs_loading = sqrt(PC1^2 + PC2^2)) %>%  
    arrange(desc(abs_loading)) %>%
    slice_max(order_by = abs_loading, n = loadings_number) %>%
    mutate(order_color = order_colors[order])
  
  # Create ggplot
  p <- ggplot() +
    geom_point(data = scores, 
               aes(x = PC1, y = PC2, 
                   fill = .data[[samples_color_metadata]],
                   shape = .data[[samples_shape_metadata]]), 
               size = 2, alpha = 0.8,
               color = "black", stroke = 0.3) +
    scale_fill_manual(values = samples_color_value) + 
    scale_shape_manual(values = c(21, 24, 23, 22, 25)) +
    new_scale_color() + 
    geom_segment(data = loadings, 
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, color = .data[[loadings_color_metadata]]),
                 arrow = arrow(length = unit(0.2, "cm")), 
                 size = 0.7, alpha = 0.9) +
    scale_color_manual(name = "Classification", values = loadings_color_value) +
    geom_text_repel(data = loadings, 
                    aes(x = PC1, y = PC2, label = .data[[loadings_taxon_level]]),
                    color = "black", size = 2, vjust = -0.5, alpha=0.7, max.overlaps = 20) +
    labs(title = "PCA Ordination Plot",
         x = pc1_label,
         y = pc2_label) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    scale_y_continuous(limits = c(-y_limit, y_limit)) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    geom_vline(xintercept = 0, color = "darkgrey") +
    theme_minimal() +
    custom_ggplot_theme +
    guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) 
  
  return(p)
}


### fit_and_analyze_model()
fit_and_analyze_model <- function(model = c("lm", "glm"),
                                  distribution = NULL,
                                  response_var,
                                  explanatory_var,
                                  data) {
  model <- match.arg(model) # Restrict to "lm" or "glm"
  # Construct the model formula
  model_formula <- as.formula(paste(response_var, "~", explanatory_var))

  # Initialize all possible return objects
  anova_result <- NULL
  md <- NULL
  simResids <- NULL

  if (model == "lm") {
    # Continuous floats (any real number)
    md <- lm(model_formula, data = data)
    anova_result <- broom::tidy(Anova(md, test.statistic = "F"))
  } else if (model == "glm") {
    # Counts (integers ≥0): poisson, quasipoisson (accounts overdispersion)
    # Proportions (0–1 continuous):	quasibinomial GLM (accounts overdispersion) but beta regression (not GLM)	Preferred over quasibinomial GLM
    # Continuous floats (any real number): gaussian, Gamma (continuous positive data, skewed to the right)
    if (is.null(distribution)) {
      stop("You must specify a distribution family for glm.")
    }

    md <- glm(model_formula, family = distribution, data = data)

    if (!grepl("^quasi", distribution)) { # Use DHARMa only on supported distributions
      simResids <- simulateResiduals(md)
      anova_result <- broom::tidy(Anova(md, test.statistic = "Wald")) # ANOVA with Chi-squared test
    } else {
      anova_result <- broom::tidy(Anova(md, test.statistic = "F")) # ANOVA with F test
    }
  }
  result <- list(
    model_fit = md,
    anova = anova_result
  )

  if (!is.null(simResids)) {
    result$simResidual <- simResids
  }

  return(result)
}

### pivot_phylo()
pivot_phylo <- function(phyloseq_obj, glom = TRUE, tax_transform = TRUE, taxon_level, tr_method) {
  if (glom == TRUE && !is.null(taxon_level)) {
    phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj)
    phyloseq_obj <- tax_glom(phyloseq_obj, taxon_level)
  } else {
    .
  }
  if (tax_transform == TRUE && !is.null(tr_method)) {
    phyloseq_obj <- tax_transform(phyloseq_obj, tr_method)
  } else {
    .
  }

  pivot_dataframe <- data.frame(otu_table(phyloseq_obj)) %>%
    rownames_to_column(var = "genome") %>%
    pivot_longer(-genome, names_to = "microsample", values_to = "abundance") %>%
    filter(abundance > 0) %>%
    left_join(data.frame(tax_table(phyloseq_obj)) %>%
      rownames_to_column(var = "genome"), by = "genome") %>%
    left_join(data.frame(sample_data(phyloseq_obj)) %>%
      rownames_to_column(var = "microsample"), by = "microsample")

  # Re-order levels
  taxa_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  # Iterate over taxonomic levels
  for (taxa in taxa_levels) {
    # Check if the column has more than one unique value
    if (taxa %in% colnames(pivot_dataframe) && length(unique(pivot_dataframe[[taxa]])) > 1) {
      # Convert each taxonomic level to a factor with levels ordered by abundance
      pivot_dataframe <- pivot_dataframe %>%
        mutate(
          !!taxa := factor(
            !!sym(taxa),
            levels = pivot_dataframe %>%
              group_by(!!sym(taxa)) %>%
              summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
              arrange(desc(total_abundance)) %>%
              pull(!!sym(taxa))
          ) # Extract ordered levels
        )
    }
  }
  return(pivot_dataframe)
}


### spatial_cryosections()
spatial_cryosections <- function(cryosection_list, metadata_df, comm_clr) {
  cryosection_dfs <- list()
  mantel_results <- list()
  mantelcor_results <- list()
  decay_dfs <- list()
  distance_decay_plots <- list()
  structure_results <- list()

  for (cryosection in cryosection_list) {
    # Filter metadata for this section
    metadata_data <- metadata_df %>%
      filter(cryosection == !!cryosection, !is.na(.data$Xcoord), !is.na(.data$Ycoord))

    # Filter community data
    comm_data <- comm_clr %>%
      data.frame() %>%
      rownames_to_column(var = "microsample") %>%
      filter(microsample %in% metadata_data$microsample) %>%
      column_to_rownames(var = "microsample")

    cryosection_dfs[[cryosection]] <- list(
      comm_clr = comm_data,
      metadata = metadata_data
    )

    # Mantel correlogram
    mantel <- vegan::mantel(
      dist(comm_data),
      dist(metadata_data[, c("Xcoord", "Ycoord")]),
      permutations = 999
    )
    mantel_results[[cryosection]] <- mantel

    # Mantel correlogram
    correlog <- vegan::mantel.correlog(
      D.eco = dist(comm_data),
      D.geo = dist(metadata_data[, c("Xcoord", "Ycoord")]),
      nperm = 999
    )
    mantelcor_results[[cryosection]] <- correlog

    # Distance decay
    toplot <- data.frame(
      spat_dist = as.numeric(dist(metadata_data[, c("Xcoord", "Ycoord")])),
      comm_dist = as.numeric(dist(comm_data))
    )
    decay_dfs[[cryosection]] <- toplot

    # Plot
    p <- ggplot(toplot, aes(x = spat_dist, y = comm_dist)) +
      geom_smooth() +
      xlab("Spatial distance (μm)") +
      ylab("Aitchison \ndistance") +
      custom_ggplot_theme +
      ggtitle(paste(cryosection))
    distance_decay_plots[[cryosection]] <- p
  }
  return(list(
    cryosection_dfs = cryosection_dfs,
    mantel_results = mantel_results,
    mantelcor_results = mantelcor_results,
    decay_dfs = decay_dfs,
    distance_decay_plots = distance_decay_plots,
    structure_results = structure_results
  ))
}

### lawsonibacter_mantel_analysis()
lawsonibacter_mantel_analysis <- function(data, circul_selection) {
  # 1. Filter and prepare presence-absence matrix
  if (circul_selection == "Y") {
    filtered_data <- data %>%
      filter(
        circul == circul_selection,
      )
  } else {
    message("circul_selection is not 'Y'; skipping circul filter.")
    filtered_data <- data
  }

  # Presence Absence
  # 1. Build presence-absence matrix
  lawsonibacter_pa <- filtered_data %>%
    mutate(presence = ifelse(abundance > 0, 1, 0)) %>%
    select(microsample, genome, presence) %>%
    pivot_wider(names_from = genome, values_from = presence, values_fill = 0) %>%
    column_to_rownames("microsample")
  # 2. Community distance
  comm_dist_pa <- vegan::vegdist(lawsonibacter_pa, method = "jaccard")
  # 3. Spatial coordinates
  coords_pa <- filtered_data %>%
    distinct(microsample, Xcoord, Ycoord) %>%
    arrange(microsample)
  coords_pa <- coords_pa[match(rownames(lawsonibacter_pa), coords_pa$microsample), ]
  # 4. Spatial distance
  spatial_dist_pa <- dist(coords_pa[, c("Xcoord", "Ycoord")])
  # 5. Mantel test
  mantel_result_pa <- vegan::mantel(comm_dist_pa, spatial_dist_pa, permutations = 999)
  # 6. Mantel correlogram
  mantelcor_result_pa <- vegan::mantel.correlog(
    D.eco = comm_dist_pa,
    D.geo = spatial_dist_pa,
    nperm = 999
  )

  # CLR
  # 1. Build count matrix and CLR transformation
  lawsonibacter_clr <- filtered_data %>%
    select(microsample, genome, abundance) %>%
    pivot_wider(names_from = genome, values_from = abundance, values_fill = 0) %>%
    column_to_rownames("microsample")
  lawsonibacter_clr <- cmultRepl(lawsonibacter_clr, method = "GBM", output = "prop", z.warning = 0.95)
  clr_transform <- function(x) {
    log(x) - mean(log(x), na.rm = TRUE)
  }
  lawsonibacter_clr <- data.frame(t(apply(lawsonibacter_clr, 1, clr_transform)))
  # 2. Community distance
  comm_dist_clr <- vegan::vegdist(lawsonibacter_clr, method = "euclidean")
  # 2. Spatial coordinates
  coords_clr <- filtered_data %>%
    distinct(microsample, Xcoord, Ycoord) %>%
    arrange(microsample) %>%
    filter(microsample %in% c(rownames(lawsonibacter_clr)))
  coords_clr <- coords_clr[match(rownames(lawsonibacter_clr), coords_clr$microsample), ]
  # 4. Spatial distance
  spatial_dist_clr <- dist(coords_clr[, c("Xcoord", "Ycoord")])
  # 5. Mantel test
  mantel_result_clr <- vegan::mantel(comm_dist_clr, spatial_dist_clr, permutations = 999)
  # 6. Mantel correlogram
  mantelcor_result_clr <- vegan::mantel.correlog(
    D.eco = comm_dist_clr,
    D.geo = spatial_dist_clr,
    nperm = 999
  )

  distance_clr_df <- data.frame(
    spat_dist = as.numeric(spatial_dist_clr),
    comm_dist = as.numeric(comm_dist_clr)
  )

  clr_lm <- aovperm(lmperm(comm_dist ~ spat_dist, data = distance_clr_df, np = 10000))

  return(list(
    # cryosection = cryosection_id,
    mantel_pa = mantel_result_pa,
    correlogram_pa = mantelcor_result_pa,
    mantel_clr = mantel_result_clr,
    correlogram_clr = mantelcor_result_clr,
    clr_lm = clr_lm
  ))
}


### ani_spatial_analysis()
ani_spatial_analysis <- function(bin_name, ani_list, ids, meta_data) {
  # Subset ANI matrix to samples of interest
  popani <- ani_list[[bin_name]][rownames(ani_list[[bin_name]]) %in% ids, ]
  popani <- popani[, colnames(popani) %in% ids] # keep only matching columns too
  diag(popani) <- NA
  popani_dist <- as.dist(1 - popani)

  # Spatial distance for those samples
  spatial_dist <- dist(meta_data %>%
    filter(microsample %in% colnames(popani)) %>%
    select(Xcoord, Ycoord))

  # Data frame for plotting and stats
  toplot <- data.frame(
    spat_dist = as.numeric(spatial_dist),
    comm_dist = as.numeric(popani_dist)
  )

  # Plot smooth relationship
  p <- ggplot(toplot, aes(x = spat_dist, y = comm_dist)) +
    geom_smooth() +
    ggtitle(paste("Spatial vs ANI distance:", bin_name))

  # print(p)

  # Permutation ANOVA test for linear relationship
  perm_aov <- aovperm(comm_dist ~ spat_dist, data = toplot, np = 10000)

  # Mantel correlogram
  mantel_res <- vegan::mantel.correlog(popani_dist, spatial_dist, nperm = 999, cutoff = TRUE, n.class = 22, r.type = "spearman")
  plot(mantel_res)
  title(main = paste("Mantel correlogram:", bin_name)) # add main title

  list(
    plot = p,
    perm_aov = perm_aov,
    mantel = mantel_res,
    spatial_dist = spatial_dist,
    mantel_summary = summary(mantel_res)
  )
}
