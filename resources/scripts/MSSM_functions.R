# format_median_iqr ()
format_median_iqr <- function(x) {
  med <- median(x, na.rm = TRUE)
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  sprintf("%.1f (%.1f–%.1f)", med, q[1], q[2])
}

# print_median_iqr()
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
  
  model <- match.arg(model)
  
  form <- stats::as.formula(
    paste(response_var, "~", explanatory_var)
  )
  
  distribution <- if (!is.null(distribution) && !is.na(distribution)) {
    tolower(trimws(as.character(distribution)))
  } else {
    NA_character_
  }
  
  if (model == "lm") {
    
    md <- stats::lm(form, data = data)
    anova_result <- car::Anova(md, type = 2, test.statistic = "F")
    test_statistic <- "F"
    
  } else {
    
    if (is.na(distribution)) {
      stop("You must specify a distribution for glm.")
    }
    
    fam <- switch(
      distribution,
      "poisson" = stats::poisson(),
      "binomial" = stats::binomial(),
      "gaussian" = stats::gaussian(),
      "gamma" = stats::Gamma(),
      "quasipoisson" = stats::quasipoisson(),
      "quasibinomial" = stats::quasibinomial(),
      stop("Unsupported GLM distribution: ", distribution)
    )
    
    md <- stats::glm(form, family = fam, data = data)
    
    if (distribution %in% c("quasipoisson", "quasibinomial")) {
      anova_result <- car::Anova(md, type = 2, test.statistic = "F")
      test_statistic <- "F"
    } else {
      anova_result <- car::Anova(md, type = 2, test.statistic = "LR")
      test_statistic <- "Chi2"
    }
  }
  
  list(
    model_fit = md,
    anova = anova_result,
    test_statistic = test_statistic,
    distribution = distribution
  )
}


### format_anova()
format_anova <- function(anova_df,
                         model_obj,
                         metric_name,
                         term_name,
                         test_statistic) {
  
  if (is.null(test_statistic) || length(test_statistic) == 0) {
    stop("test_statistic is missing. Check run_metric() is passing fit$test_statistic.")
  }
  
  format_p_value <- function(p) {
    ifelse(
      is.na(p),
      NA_character_,
      ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
    )
  }
  
  df <- as.data.frame(anova_df)
  
  if (!"term" %in% names(df)) {
    df$term <- rownames(df)
  }
  
  row <- df[df$term == term_name, , drop = FALSE]
  
  if (nrow(row) == 0) {
    return(NA_character_)
  }
  
  df_col <- if ("df" %in% names(row)) {
    "df"
  } else {
    "Df"
  }
  
  if (test_statistic == "F") {
    
    stat_col <- if ("statistic" %in% names(row)) {
      "statistic"
    } else {
      grep("^F value$|^F$", names(row), value = TRUE)[1]
    }
    
    p_col <- if ("p.value" %in% names(row)) {
      "p.value"
    } else {
      grep("Pr\\(>F\\)", names(row), value = TRUE)[1]
    }
    
    if (is.na(stat_col) || is.na(p_col) || is.na(df_col)) {
      stop("Could not find F statistic, p-value, or df column. Columns are: ",
           paste(names(row), collapse = ", "))
    }
    
    sprintf(
      "%s: F(%s, %s) = %.2f, p = %s for %s",
      metric_name,
      row[[df_col]],
      stats::df.residual(model_obj),
      row[[stat_col]],
      format_p_value(row[[p_col]]),
      term_name
    )
    
  } else {
    
    stat_col <- if ("statistic" %in% names(row)) {
      "statistic"
    } else {
      grep("Chisq|LR", names(row), value = TRUE)[1]
    }
    
    p_col <- if ("p.value" %in% names(row)) {
      "p.value"
    } else {
      grep("Pr\\(>Chisq\\)", names(row), value = TRUE)[1]
    }
    
    if (is.na(stat_col) || is.na(p_col) || is.na(df_col)) {
      stop("Could not find Chi-square statistic, p-value, or df column. Columns are: ",
           paste(names(row), collapse = ", "))
    }
    
    sprintf(
      "%s: χ²(%s) = %.2f, p = %s for %s",
      metric_name,
      row[[df_col]],
      row[[stat_col]],
      format_p_value(row[[p_col]]),
      term_name
    )
  }
}


### run_metric()
run_metric <- function(metric,
                       response_var,
                       model,
                       distribution = NULL,
                       section_name = NULL,
                       explanatory_var,
                       term_name,
                       data) {
  
  if (!is.null(section_name)) {
    data <- data %>%
      dplyr::filter(section == section_name)
  }
  
  fit <- fit_and_analyze_model(
    model = model,
    distribution = distribution,
    response_var = response_var,
    explanatory_var = explanatory_var,
    data = data
  )
  
  tibble::tibble(
    metric = metric,
    section_name = section_name,
    report = format_anova(
      anova_df = fit$anova,
      model_obj = fit$model_fit,
      metric_name = metric,
      term_name = term_name,
      test_statistic = fit$test_statistic
    )
  )
}
###

### pivot_phylo()
pivot_phylo <- function(phyloseq_obj, glom = TRUE, tax_transform = TRUE, taxon_level, tr_method) {
  if (glom == TRUE && !is.null(taxon_level)) {
    phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj)
    phyloseq_obj <- tax_glom(phyloseq_obj, taxon_level)
  } else {
    .
  }
  if (tax_transform == TRUE && !is.null(tr_method)) {
    phyloseq_obj <- microbiome::transform(phyloseq_obj, tr_method)
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
  
  taxa_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  for (taxa in taxa_levels) {
    if (taxa %in% colnames(pivot_dataframe) && length(unique(pivot_dataframe[[taxa]])) > 1) {
      pivot_dataframe <- pivot_dataframe %>%
        mutate(
          !!taxa := factor(
            !!sym(taxa),
            levels = pivot_dataframe %>%
              group_by(!!sym(taxa)) %>%
              summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
              arrange(desc(total_abundance)) %>%
              pull(!!sym(taxa))
          ) 
        )
    }
  }
  return(pivot_dataframe)
}


### compare_taxa_overlap()
compare_taxa_overlap <- function(data,
                                 group_var,
                                 group_a,
                                 group_b,
                                 facet_var = NULL,
                                 taxon_var = "species",
                                 abundance_var = "abundance",
                                 sample_var = "microsample") {
  
  df <- data %>%
    filter(.data[[group_var]] %in% c(group_a, group_b)) %>%
    filter(!is.na(.data[[taxon_var]])) %>%
    mutate(
      group = .data[[group_var]],
      taxon = .data[[taxon_var]],
      sample = .data[[sample_var]],
      abundance = .data[[abundance_var]],
      present = abundance > 0
    )
  
  ## Check facet variable exists
  if (!is.null(facet_var) && !facet_var %in% names(df)) {
    stop("'", facet_var, "' is not a column in 'data'.")
  }
  
  ## Detection (facet-aware if facet_var is provided)
  detection_vars <- c(
    if (!is.null(facet_var)) facet_var,
    "group",
    "taxon"
  )
  
  detection <- df %>%
    filter(present) %>%
    distinct(across(all_of(detection_vars))) %>%
    mutate(present = TRUE) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = present,
      values_fill = FALSE
    ) %>%
    mutate(
      Detection = case_when(
        .data[[group_a]] & !.data[[group_b]] ~ paste(group_a, "only"),
        !.data[[group_a]] & .data[[group_b]] ~ paste(group_b, "only"),
        .data[[group_a]] & .data[[group_b]] ~ "Shared",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Detection))
  
  ## Number of samples per group (and facet)
  sample_total_vars <- c(
    if (!is.null(facet_var)) facet_var,
    "group",
    "sample"
  )
  
  sample_count_vars <- c(
    if (!is.null(facet_var)) facet_var,
    "group"
  )
  
  sample_totals <- df %>%
    distinct(across(all_of(sample_total_vars))) %>%
    group_by(across(all_of(sample_count_vars))) %>%
    summarise(
      n_samples = n(),
      .groups = "drop"
    )
  
  ## Mean abundance and prevalence
  plot_group_vars <- c(
    if (!is.null(facet_var)) facet_var,
    "group",
    "taxon"
  )
  
  plot_df <- df %>%
    group_by(across(all_of(plot_group_vars))) %>%
    summarise(
      mean_abundance = mean(abundance, na.rm = TRUE),
      prevalence_raw = sum(present),
      .groups = "drop"
    ) %>%
    left_join(
      sample_totals,
      by = sample_count_vars
    ) %>%
    mutate(
      prevalence = prevalence_raw / n_samples
    ) %>%
    left_join(
      detection,
      by = c(
        if (!is.null(facet_var)) facet_var,
        "taxon"
      )
    ) %>%
    filter(!is.na(Detection))
  
  list(
    plot_df = plot_df,
    detection = detection,
    sample_totals = sample_totals
  )
}

### prepare_spatial_section()
prepare_spatial_section <- function(cryosection_id,
                                    comm_data,
                                    metadata,
                                    z.warning = 0.95) {
  
  comm_section <- comm_data %>%
    as.data.frame() %>%
    rownames_to_column(var = "microsample") %>%
    left_join(metadata, by = "microsample") %>%
    filter(cryosection == cryosection_id) %>%
    column_to_rownames(var = "microsample") %>%
    select(contains("bin_"))
  
  comm_zerRepl <- cmultRepl(
    comm_section,
    method = "GBM",
    output = "prop",
    z.warning = z.warning
  )
  
  metadata_section <- metadata %>%
    filter(microsample %in% rownames(comm_zerRepl)) %>%
    mutate(
      Xcoord = as.numeric(Xcoord),
      Ycoord = as.numeric(Ycoord)
    )
  
  comm_red <- comm_section %>%
    select(all_of(colnames(comm_zerRepl))) %>%
    filter(rownames(comm_section) %in% rownames(comm_zerRepl))
  
  list(
    comm = comm_section,
    comm_zerRepl = comm_zerRepl,
    metadata_section = metadata_section,
    comm_red = comm_red
  )
}

