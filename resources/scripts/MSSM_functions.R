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
  desaturate_palette <- function(cols, amount = 0.12) {
    desaturated <- vapply(cols, function(col) {
      rgb_col <- grDevices::col2rgb(col) / 255
      grey_val <- mean(rgb_col)
      mixed <- rgb_col * (1 - amount) + grey_val * amount
      grDevices::rgb(mixed[1], mixed[2], mixed[3])
    }, character(1))
    stats::setNames(desaturated, names(cols))
  }
  make_breaks <- function(x, n = 5) {
    rng <- range(x, na.rm = TRUE)
    brks <- pretty(rng, n = n)
    brks <- brks[brks >= rng[1] & brks <= rng[2]]
    if (length(brks) < 2) brks <- unique(rng)
    brks
  }
  genome_size_colors <- c(
    "#eef3f6",
    "#c9d7e1",
    "#98b2c2",
    "#66879c",
    "#35556b"
  )
  contamination_colors <- c(
    "#f6f1ed",
    "#e2cfc3",
    "#c7a08d",
    "#a36d58",
    "#6f3f32"
  )
  completeness_colors <- c(
    "#eef2ea",
    "#cfd9c6",
    "#a8bb9c",
    "#6f8f6e",
    "#3f5f45"
  )
  genome_metadata_ordered <- genome_metadata %>%
    arrange(match(genome, tree$tip.label)) %>%
    mutate(
      phylum = factor(phylum, levels = unique(phylum)),
      order = factor(order, levels = unique(order))
    )
  circularised_genomes <- genome_metadata_ordered %>%
    filter(as.character(circularity) %in% c("1", "TRUE", "true")) %>%
    pull(genome)
  completeness_breaks <- make_breaks(genome_metadata_ordered$completeness, n = 5)
  contamination_breaks <- make_breaks(genome_metadata_ordered$contamination, n = 5)
  genome_size_breaks <- make_breaks(genome_metadata_ordered$length, n = 5)
  
  ### basal tree
  circular_tree_base <- force.ultrametric(tree, method = "extend") # extend to ultrametric for visualisation
  circular_tree_base$edge.length <- circular_tree_base$edge.length * 2.05
  circular_tree <- ggtree(circular_tree_base, layout = "fan", open.angle = 10, size = 0.35, color = "grey55")
  circularised_tip_data <- circular_tree$data %>%
    filter(isTip, label %in% circularised_genomes)
  circular_tree <- circular_tree +
    geom_tippoint(data = circularised_tip_data, color = "#fcbb6d", size = 1.2, alpha = 0.95)
  
  ### completeness ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_stepsn(
      colors = completeness_colors,
      breaks = completeness_breaks,
      name = "Completeness %",
      na.value = "white",
      guide = guide_colorsteps(
        barwidth = grid::unit(4.2, "cm"),
        barheight = grid::unit(0.35, "cm")
      )
    ) +
    geom_fruit(data=genome_metadata_ordered, geom=geom_bar, 
               mapping = aes(x=completeness, y=genome, fill=completeness),
               offset = 0.05, pwidth = 0.07, orientation="y", stat="identity")
  
  ### contamination ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_stepsn(
      colors = contamination_colors,
      breaks = contamination_breaks,
      name = "Contamination %",
      na.value = "white",
      guide = guide_colorsteps(
        barwidth = grid::unit(4.2, "cm"),
        barheight = grid::unit(0.35, "cm")
      )
    ) +
    geom_fruit(data=genome_metadata_ordered, geom=geom_bar,
               mapping = aes(x=contamination, y=genome, fill=contamination),
               offset = 0.015, pwidth = 0.07, orientation="y", stat="identity")
  
  ### genome size ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_stepsn(
      colors = genome_size_colors,
      breaks = genome_size_breaks,
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      name = "Genome size",
      na.value = "white",
      guide = guide_colorsteps(
        barwidth = grid::unit(4.2, "cm"),
        barheight = grid::unit(0.35, "cm")
      )
    ) +
    geom_fruit(data = genome_metadata_ordered, geom = geom_bar,
               mapping = aes(x = length, y = genome, fill = length),
               offset = 0.015, pwidth = 0.07, orientation = "y", stat = "identity")
  
  circular_tree <- circular_tree +
    theme(plot.margin = margin(0, 0, 0, 0), 
          panel.margin = margin(0, 0, 0, 0))
  
  ### order ring 
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_manual(
      values = desaturate_palette(order_colors, amount = 0.12),
      name = "Order",
      guide = guide_legend(ncol = 2, byrow = FALSE)
    ) +
    geom_fruit(data = genome_metadata_ordered, geom = geom_bar,
               mapping = aes(x = 1, y = genome, fill = order),
               offset = 0.015, pwidth = 0.048, orientation = "y", stat = "identity",
               width = 1.02, linewidth = 0, color = NA)
  
  ### phylum ring
  circular_tree <- circular_tree +
    new_scale_fill() +
    scale_fill_manual(values = desaturate_palette(phylum_colors, amount = 0.12), name = "Phylum") +
    geom_fruit(data = genome_metadata_ordered, geom = geom_bar,
               mapping = aes(x = 1, y = genome, fill = phylum),
               offset = 0.001, pwidth = 0.048, orientation = "y", stat = "identity",
               width = 1.03, linewidth = 0, color = NA) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.05, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.key.width = grid::unit(0.5, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # text
  circular_tree <-  circular_tree +
    annotate('text', x=4.00, y=0, label='Circularised', size=3.5, hjust = 0) +
    annotate('text', x=4.40, y=0, label='Completeness', size=3.5, hjust = 0) +
    annotate('text', x=4.80, y=0, label='Contamination', size=3.5, hjust = 0) +
    annotate('text', x=5.10, y=0, label='Genome size', size=3.5, hjust = 0) +
    annotate('text', x=5.45, y=0, label='Order', size=3.5, hjust = 0) +
    annotate('text', x=5.70, y=0, label='Phylum', size=3.5, hjust = 0) +
    coord_cartesian(clip = "off")
  
  # circular tree
  p <-circular_tree %>% open_tree(55) %>% rotate_tree(90) 
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
plot_pca <- function(df, 
                     samples_color_metadata, samples_shape_metadata, 
                     samples_color_value, loadings_color_metadata, 
                     loadings_color_value, loadings_taxon_level,
                     sample_metadata, genome_metadata, order_colors, 
                     custom_ggplot_theme, scaling_factor_value = 1.5, 
                     loadings_number = 10,
                     group_overlay = c("none", "ellipse", "kde", "both"),
                     overlay_group_metadata = NULL,
                     overlay_color_values = NULL,
                     overlay_legend = FALSE,
                     overlay_legend_title = "Group overlay",
                     ellipse_level = 0.95,
                     ellipse_linewidth = 0.55,
                     ellipse_alpha = 0.9,
                     kde_contour_type = c("line", "filled"),
                     kde_bins = 4,
                     kde_breaks = NULL,
                     kde_contour_var = c("density", "ndensity"),
                     kde_adjust = c(1, 1),
                     kde_n = 100,
                     kde_linewidth = 0.6,
                     kde_alpha = 0.9,
                     plot_outliers = TRUE,
                     kde_exclude_outliers_from_fit = FALSE,
                     percent_outliers = 0.1,
                     outlier_shape = 4,
                     outlier_size = 1.8) {
  
  group_overlay <- match.arg(group_overlay)
  kde_contour_type <- match.arg(kde_contour_type)
  kde_contour_var <- match.arg(kde_contour_var)
  percent_outliers <- max(0, min(0.49, percent_outliers))
  
  scores <- rownames_to_column(as.data.frame(df$x), var = "microsample")
  scores <- left_join(scores, sample_metadata, by = join_by(microsample == microsample))
  
  x_limit <- max(abs(scores$PC1))
  y_limit <- max(abs(scores$PC2))
  
  variance_explained <- (df$sdev^2) / sum(df$sdev^2) * 100
  pc1_label <- paste0("PC1: ", round(variance_explained[1], 2), "% variance explained")
  pc2_label <- paste0("PC2: ", round(variance_explained[2], 2), "% variance explained")
  
  scaling_factor <- scaling_factor_value
  
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
  
  shape_levels <- unique(as.character(stats::na.omit(scores[[samples_shape_metadata]])))
  shape_values <- setNames(
    rep(c(21, 24, 23, 22, 25), length.out = length(shape_levels)),
    shape_levels
  )
  
  p <- ggplot() +
    geom_point(data = scores, 
               aes(x = PC1, y = PC2, 
                   fill = .data[[samples_color_metadata]],
                   shape = .data[[samples_shape_metadata]]), 
               size = 2, alpha = 0.8,
               color = "black", stroke = 0.3) +
    scale_fill_manual(
      values = samples_color_value,
      name = samples_color_metadata,
      drop = FALSE,
      guide = guide_legend(override.aes = list(shape = 21, color = "black"))
    ) + 
    scale_shape_manual(values = shape_values, drop = FALSE)
  
  use_overlay <- group_overlay != "none"
  overlay_data <- NULL
  overlay_kde_data <- NULL
  overlay_colors <- NULL
  overlay_outliers <- tibble::tibble()
  
  if (use_overlay) {
    if (is.null(overlay_group_metadata) || !overlay_group_metadata %in% colnames(scores)) {
      warning("Overlay skipped: 'overlay_group_metadata' is missing in sample metadata.")
      use_overlay <- FALSE
    } else {
      overlay_data <- scores %>%
        filter(!is.na(PC1), !is.na(PC2), !is.na(.data[[overlay_group_metadata]])) %>%
        mutate(overlay_group = factor(.data[[overlay_group_metadata]]))
      
      if (nrow(overlay_data) < 4 || nlevels(overlay_data$overlay_group) < 2) {
        warning("Overlay skipped: insufficient points/groups for contour/ellipse overlay.")
        use_overlay <- FALSE
      } else {
        overlay_levels <- levels(overlay_data$overlay_group)
        fallback_cols <- setNames(grDevices::hcl.colors(length(overlay_levels), palette = "Dark 3"), overlay_levels)
        
        if (!is.null(overlay_color_values)) {
          if (is.null(names(overlay_color_values))) {
            overlay_colors <- setNames(rep(overlay_color_values, length.out = length(overlay_levels)), overlay_levels)
          } else {
            overlay_colors <- overlay_color_values[overlay_levels]
            names(overlay_colors) <- overlay_levels
          }
        } else if (!is.null(samples_color_value)) {
          if (is.null(names(samples_color_value))) {
            overlay_colors <- setNames(rep(samples_color_value, length.out = length(overlay_levels)), overlay_levels)
          } else {
            overlay_colors <- samples_color_value[overlay_levels]
            names(overlay_colors) <- overlay_levels
          }
        }
        
        if (is.null(overlay_colors)) {
          overlay_colors <- fallback_cols
        } else {
          na_idx <- is.na(overlay_colors)
          overlay_colors[na_idx] <- fallback_cols[names(overlay_colors)[na_idx]]
        }
        
        if (group_overlay %in% c("kde", "both") && (plot_outliers || kde_exclude_outliers_from_fit)) {
          overlay_outliers <- overlay_data %>%
            group_by(overlay_group) %>%
            group_modify(~{
              if (nrow(.x) < 5) {
                return(.x[0, , drop = FALSE])
              }
              cov_mat <- tryCatch(stats::cov(.x[, c("PC1", "PC2")]), error = function(e) NA)
              if (anyNA(cov_mat) || abs(det(cov_mat)) < .Machine$double.eps) {
                return(.x[0, , drop = FALSE])
              }
              d2 <- stats::mahalanobis(.x[, c("PC1", "PC2")], center = colMeans(.x[, c("PC1", "PC2")]), cov = cov_mat)
              cutoff <- stats::quantile(d2, probs = 1 - percent_outliers, na.rm = TRUE)
              .x[d2 > cutoff, , drop = FALSE]
            }) %>%
            ungroup()
          
          if (kde_exclude_outliers_from_fit && nrow(overlay_outliers) > 0) {
            overlay_kde_data <- overlay_data %>%
              anti_join(
                overlay_outliers %>% select(microsample, overlay_group),
                by = c("microsample", "overlay_group")
              )
          }
        }
      }
    }
  }
  
  if (use_overlay) {
    has_ellipse <- group_overlay %in% c("ellipse", "both")
    has_kde <- group_overlay %in% c("kde", "both")
    use_overlay_color_scale <- FALSE
    if (is.null(overlay_kde_data)) overlay_kde_data <- overlay_data
    
    if (has_ellipse) {
      p <- p + stat_ellipse(
        data = overlay_data,
        aes(x = PC1, y = PC2, color = overlay_group, group = overlay_group),
        type = "norm",
        level = ellipse_level,
        linewidth = ellipse_linewidth,
        alpha = ellipse_alpha,
        show.legend = overlay_legend
      )
      use_overlay_color_scale <- TRUE
    }
    
    if (has_kde) {
      if (kde_contour_type == "line") {
        p <- p + stat_density_2d(
          data = overlay_kde_data,
          aes(x = PC1, y = PC2, color = overlay_group, group = overlay_group),
          contour = TRUE,
          contour_var = kde_contour_var,
          bins = kde_bins,
          breaks = kde_breaks,
          adjust = kde_adjust,
          n = kde_n,
          linewidth = kde_linewidth,
          alpha = kde_alpha,
          show.legend = overlay_legend
        )
        use_overlay_color_scale <- TRUE
      } else {
        p <- p +
          ggnewscale::new_scale_fill() +
          stat_density_2d(
            data = overlay_kde_data,
            aes(x = PC1, y = PC2, fill = overlay_group, group = overlay_group, alpha = after_stat(level)),
            geom = "polygon",
            contour = TRUE,
            contour_var = kde_contour_var,
            bins = kde_bins,
            breaks = kde_breaks,
            adjust = kde_adjust,
            n = kde_n,
            color = NA,
            show.legend = overlay_legend && !has_ellipse
          ) +
          scale_fill_manual(
            values = overlay_colors,
            name = overlay_legend_title,
            guide = if (overlay_legend && !has_ellipse) "legend" else "none"
          ) +
          scale_alpha_continuous(range = c(0.10, kde_alpha), guide = "none") +
          stat_density_2d(
            data = overlay_kde_data,
            aes(x = PC1, y = PC2, color = overlay_group, group = overlay_group),
            contour = TRUE,
            contour_var = kde_contour_var,
            bins = kde_bins,
            breaks = kde_breaks,
            adjust = kde_adjust,
            n = kde_n,
            linewidth = kde_linewidth * 0.7,
            alpha = pmin(1, kde_alpha + 0.05),
            show.legend = FALSE
          )
      }
      
      if (plot_outliers && nrow(overlay_outliers) > 0) {
        if (kde_contour_type == "line" || has_ellipse) {
          p <- p + geom_point(
            data = overlay_outliers,
            aes(x = PC1, y = PC2, color = overlay_group),
            shape = outlier_shape,
            size = outlier_size,
            stroke = 0.5,
            alpha = 0.9,
            show.legend = FALSE
          )
          use_overlay_color_scale <- TRUE
        } else {
          p <- p + geom_point(
            data = overlay_outliers,
            aes(x = PC1, y = PC2),
            shape = outlier_shape,
            size = outlier_size,
            stroke = 0.5,
            alpha = 0.9,
            color = "black",
            show.legend = FALSE
          )
        }
      }
    }
    
    if (use_overlay_color_scale) {
      p <- p + scale_color_manual(values = overlay_colors, name = overlay_legend_title)
    }
  }
  
  if (use_overlay) {
    p <- p +
      ggnewscale::new_scale_fill() +
      scale_fill_manual(values = samples_color_value, guide = "none", drop = FALSE) +
      geom_point(
        data = scores,
        aes(x = PC1, y = PC2,
            fill = .data[[samples_color_metadata]],
            shape = .data[[samples_shape_metadata]]),
        size = 2, alpha = 0.8,
        color = "black", stroke = 0.3,
        show.legend = FALSE
      )
  }
  
  p <- p +
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
    custom_ggplot_theme
  
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

# prepare_spatial_section()
prepare_spatial_section <- function(cryosection_id,
                                    comm_data,
                                    metadata,
                                    z.warning = 0.95) {
  
  comm <- comm_data %>%
    as.data.frame() %>%
    rownames_to_column(var = "microsample") %>%
    bind_cols(metadata) %>%
    filter(cryosection == cryosection_id) %>%
    column_to_rownames(var = "microsample...1") %>%
    select(contains("bin_"))
  
  comm_zerRepl <- cmultRepl(
    comm,
    method = "GBM",
    output = "prop",
    z.warning = z.warning
  )
  
  metadata_section <- metadata %>%
    filter(microsample %in% rownames(comm_zerRepl))
  
  comm_red <- comm %>%
    select(colnames(comm_zerRepl)) %>%
    filter(rownames(comm) %in% rownames(comm_zerRepl))
  
  list(
    comm = comm,
    comm_zerRepl = comm_zerRepl,
    metadata = metadata_section,
    comm_red = comm_red
  )
}

## 3. Spatial RLQ function
run_spatial_rlq <- function(cryosection_id,
                            comm_red,
                            metadata_section,
                            gift_pcoa,
                            genome_tree,
                            root = 0.5,
                            correlogram_order = 8) {

  comp <- decostand(comm_red, MARGIN = 1, method = "total")
  comp <- comp^root
  colnames(comp) <- gsub("\\.", ":", colnames(comp))
  
  env <- data.frame(
    log_seq_counts = log(metadata_section$after_filtering.total_bases),
    div = metadata_section$richness,
    host_dist = log(metadata_section$distance_host)
  )
  
  funct_PCOA <- gift_pcoa$vectors[
    rownames(gift_pcoa$vectors) %in% colnames(comp),
    1:2,
    drop = FALSE
  ]
  
  phy <- drop.tip(
    genome_tree,
    setdiff(genome_tree$tip.label, rownames(funct_PCOA))
  )
  
  spa <- metadata_section[, c("Xcoord", "Ycoord")]
  
  comp <- comp[, match(phy$tip.label, colnames(comp)), drop = FALSE]
  
  funct_PCOA <- data.frame(
    funct_PCOA[
      match(phy$tip.label, rownames(funct_PCOA)),
      ,
      drop = FALSE
    ]
  )
  
  alignment_check <- c(
    phy_comp = mean(phy$tip.label == colnames(comp)),
    phy_funct = mean(phy$tip.label == rownames(funct_PCOA))
  )
  
  phylog <- newick2phylog(write.tree(phy))
  
  colnames(comp) <- gsub(":", "_", colnames(comp))
  rownames(funct_PCOA) <- gsub(":", "_", rownames(funct_PCOA))
  
  coacomp <- dudi.coa(comp, scan = FALSE, nf = ncol(comp))
  
  nb1 <- graph2nb(gabrielneigh(as.matrix(spa)), sym = TRUE)
  
  correlograms <- list(
    div = sp.correlogram(
      nb1,
      log(env$div),
      order = correlogram_order,
      method = "I"
    ),
    log_seq_counts = sp.correlogram(
      nb1,
      env$log_seq_counts,
      order = correlogram_order,
      method = "I"
    ),
    host_dist = sp.correlogram(
      nb1,
      env$host_dist,
      order = correlogram_order,
      method = "I"
    )
  )
  
  lw1 <- nb2listw(nb1)
  
  nb1_neigh <- nb2neig(nb1)
  vecspa <- scores.neig(nb1_neigh)
  
  pcaspa <- dudi.pca(
    vecspa,
    row.w = coacomp$lw,
    scan = FALSE,
    nf = ncol(vecspa)
  )
  
  pcaenv <- dudi.pca(
    env,
    row.w = coacomp$lw,
    scannf = FALSE,
    nf = 2
  )
  
  pcophy <- dudi.pco(
    as.dist(as.matrix(phylog$Wdist)[names(comp), names(comp)]),
    coacomp$cw,
    full = TRUE
  )
  
  disT <- dist.ktab(
    ktab.list.df(list(funct_PCOA)),
    c("Q"),
    scan = FALSE
  )
  
  pcotraits <- dudi.pco(
    disT,
    coacomp$cw,
    full = TRUE
  )
  
  rlqmix <- rlqESLTP(
    pcaenv,
    pcaspa,
    coacomp,
    pcotraits,
    pcophy,
    scan = FALSE,
    nf = 2
  )
  
  eig_prop <- rlqmix$eig / sum(rlqmix$eig)
  
  rlq_scores <- data.frame(
    microsample = rownames(rlqmix$lR),
    rlqmix$lR
  )
  
  rlq_scores$microsample <- rownames(comm_red)[
    as.numeric(rlq_scores$microsample)
  ]
  
  rlq_scores <- rlq_scores %>%
    left_join(
      metadata_section %>%
        select(microsample, Xcoord, Ycoord, after_filtering.total_bases,
               richness,distance_host) %>% 
        mutate(log_seq_counts= log(after_filtering.total_bases),
               log_distance_host = log(metadata_section$distance_host)),
      by = "microsample"
    )
  
  list(
    cryosection = cryosection_id,
    comp = comp,
    env = env,
    funct_PCOA = funct_PCOA,
    phy = phy,
    phylog = phylog,
    spa = spa,
    coacomp = coacomp,
    nb1 = nb1,
    correlograms = correlograms,
    lw1 = lw1,
    pcaspa = pcaspa,
    pcaenv = pcaenv,
    pcophy = pcophy,
    disT = disT,
    pcotraits = pcotraits,
    rlqmix = rlqmix,
    eig_prop = eig_prop,
    rlq_scores = rlq_scores,
    alignment_check = alignment_check
  )
}
