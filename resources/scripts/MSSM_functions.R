# load_working_data()
# The stored workspace was created with save.image(), so it also embeds plotting themes and copies of helper functions. 
# Themes serialised by a newer ggplot2 break theme merging; functions overwrite the versions sourced from MSSM_functions.R / JEC_1743_sm_apps5.R (and dplyr::select / dplyr::recode set in index.Rmd). 
# We skip both so index.Rmd remains the source of truth for code.
# Requires those scripts (and the dplyr aliases) to be sourced before this runs, which is the normal bookdown order.
load_working_data <- function(path = "resources/working_data_object.Rdata",
                              envir = globalenv()) {
  data_env <- new.env()
  load(path, envir = data_env)
  objects <- ls(data_env, all.names = TRUE)
  skip <- vapply(objects, function(x) {
    obj <- get(x, envir = data_env)
    inherits(obj, "theme") || is.function(obj)
  }, logical(1))
  list2env(mget(objects[!skip], envir = data_env), envir = envir)
  invisible(objects[skip])
}


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

## Spatial RLQ function
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

## Spatial RLQ tables

# match_genomes()
# RLQ tables name genomes with "_" (GPB_bin_000056), the tree and the GIFT tables with ":" (GPB:bin_000056).
match_genomes <- function(norm_names, pool) {
  pool[gsub(":", "_", pool) %in% norm_names]
}

# gifts_subset()
gifts_subset <- function(norm_names, gifts = macromag_gifts) {
  gifts[match_genomes(norm_names, rownames(gifts)), , drop = FALSE]
}

# tree_subset()
tree_subset <- function(norm_names, tree = macromag_tree) {
  ape::keep.tip(tree, match_genomes(norm_names, tree$tip.label))
}

# safe_rotate()
# Node numbers refer to a tree pruned to a given genome set, so rotate only when the node is an internal node of that tree.
safe_rotate <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node > n_tip && node <= n_tip + tree$Nnode) return(ape::rotate(tree, node))
  warning("node ", node, " is not an internal node of this tree", call. = FALSE)
  tree
}

# rlq_traits_table()
# Genome-side (P+T) loadings of one cryosection: one row per genome, one AxcQ<axis>_<cryosection> column per requested axis.
rlq_traits_table <- function(section, axes = 1:2,
                             results = spatial_rlq_results_flipped) {
  lq <- as.data.frame(results[[section]]$rlqmix$lQ)
  out <- tibble(genome = rownames(lq))
  for (ax in axes) {
    out[[paste0("AxcQ", ax, "_", section)]] <- lq[[paste0("AxcQ", ax)]]
  }
  out
}

# rlq_space_table()
# Sample-side (E+S) scores of one cryosection: coordinates, one AxcR<axis>_<cryosection> column per axis and the environmental RLQ inputs.
rlq_space_table <- function(section, axes = 1:2,
                            results = spatial_rlq_results_flipped,
                            datasets = spatial_datasets) {
  scores <- results[[section]]$rlq_scores
  out <- scores %>% select(microsample, Xcoord, Ycoord)
  for (ax in axes) {
    out[[paste0("AxcR", ax, "_", section)]] <- scores[[paste0("AxcR", ax)]]
  }
  out %>%
    left_join(
      datasets[[section]]$metadata %>%
        transmute(
          microsample,
          div = richness,
          log_seq_counts = log(after_filtering.total_bases),
          distance_host
        ),
      by = "microsample"
    ) %>%
    as_tibble()
}

# rlq_repl_traits_table()
# Genome-side loadings of several cryosections side by side, restricted to the genomes they share.
rlq_repl_traits_table <- function(sections, axis,
                                  results = spatial_rlq_results_flipped) {
  Reduce(
    function(x, y) inner_join(x, y, by = "genome"),
    lapply(sections, rlq_traits_table, axes = axis, results = results)
  )
}

# make_extRLQ_tree_data()
# Per-tip annotation: order/family plus the genome-side axis columns.
make_extRLQ_tree_data <- function(tree, traits,
                                  genome_metadata = macromag_genomemetadata) {
  tibble(genome = tree$tip.label) %>%
    mutate(norm = gsub(":", "_", genome)) %>%
    left_join(genome_metadata %>% select(genome, order, family), by = "genome") %>%
    left_join(traits, by = c("norm" = "genome"))
}


## Spatial RLQ figures

# extRLQ_rdbu_scale()
# Shared diverging scale of all sample-side RLQ plots.
extRLQ_rdbu_scale <- function(lim, name) {
  scale_color_fermenter(palette = "RdBu", n.breaks = 7,
                        limits = c(-lim, lim), name = name)
}

# extRLQ_axis_r_lim()
# Symmetric colour limit spanning several cryosections.
extRLQ_axis_r_lim <- function(space_dfs, axis_cols) {
  vals <- unlist(Map(function(df, col) df[[col]], space_dfs, axis_cols))
  max(abs(vals), na.rm = TRUE)
}

# build_extRLQ_grad_tree()
# Tree with tip points coloured by an RLQ axis (`axis_col`), and clade boxes and clade labels at `tax_rank` drawn as a fade of the taxon colour.
build_extRLQ_grad_tree <- function(tree, tree_data, axis_col, axis_label,
                                   tax_rank = "order") {
  color_map <- switch(tax_rank,
                      family = family_colors,
                      order  = order_colors,
                      stop("tax_rank must be 'family' or 'order'"))

  axis_lim <- max(abs(tree_data[[axis_col]]), na.rm = TRUE)

  grad <- ggtree(tree, size = 0.6, ladderize = FALSE) %<+% tree_data
  tips <- grad$data %>% filter(isTip)
  x_root <- min(grad$data$x, na.rm = TRUE)
  label_x <- max(tips$x, na.rm = TRUE) + 0.04

  labels <- tips %>%
    filter(!is.na(.data[[tax_rank]])) %>%
    mutate(taxon = .data[[tax_rank]]) %>%
    arrange(desc(y)) %>%
    distinct(taxon, .keep_all = TRUE) %>%
    transmute(taxon, x = label_x, y)

  boxes <- tips %>%
    filter(!is.na(.data[[tax_rank]])) %>%
    mutate(taxon = .data[[tax_rank]]) %>%
    group_by(taxon) %>%
    summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")

  grad_layers <- lapply(seq_len(nrow(boxes)), function(i) {
    tax <- boxes$taxon[i]
    tax_fill <- if (tax %in% names(color_map)) color_map[[tax]] else "grey80"
    tax_grad <- grid::linearGradient(
      colours = c(scales::alpha(tax_fill, 0), scales::alpha(tax_fill, 0.65)),
      x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(fill = tax_grad, col = NA)),
      xmin = x_root, xmax = label_x,
      ymin = boxes$ymin[i], ymax = boxes$ymax[i])
  })

  grad <- grad +
    grad_layers +
    geom_tippoint(aes(fill = .data[[axis_col]]),
                  shape = 21, size = 2.8, stroke = 0.4, color = "grey30") +
    scale_fill_fermenter(palette = "RdBu", n.breaks = 7,
                         limits = c(-axis_lim, axis_lim), name = axis_label) +
    geom_text2(data = labels, mapping = aes(x = x, y = y, label = taxon),
               inherit.aes = FALSE, size = 2.6, hjust = 0, vjust = 0.5) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title.position = "top",
          plot.margin = margin(5.5, 90, 5.5, 5.5, "pt"))

  # gradient boxes to the back, so branches and tips sit on top
  is_box <- vapply(grad$layers,
                   function(l) inherits(l$geom, "GeomCustomAnn"), logical(1))
  grad$layers <- c(grad$layers[is_box], grad$layers[!is_box])
  grad
}

# build_extRLQ_family_tree_bars()
# Tree with tip points coloured by family (which frees the fill aesthetic) and one bar panel per column in `bar_cols`, all on a single shared RdBu scale.
build_extRLQ_family_tree_bars <- function(tree, tree_data, bar_cols, bar_legend,
                                         tax_rank = "order",
                                         bar_pwidth = 0.25,
                                         bar_offset = 0.18,
                                         first_offset = 0.5,
                                         bar_width = 0.7) {
  color_map <- switch(tax_rank,
                      family = family_colors,
                      order  = order_colors,
                      stop("tax_rank must be 'family' or 'order'"))

  bar_lim <- max(abs(unlist(tree_data[bar_cols])), na.rm = TRUE)

  grad <- ggtree(tree, size = 0.6, ladderize = FALSE) %<+% tree_data
  tips <- grad$data %>% filter(isTip)
  x_root <- min(grad$data$x, na.rm = TRUE)
  label_x <- max(tips$x, na.rm = TRUE) + 0.04

  labels <- tips %>%
    filter(!is.na(.data[[tax_rank]])) %>%
    mutate(taxon = .data[[tax_rank]]) %>%
    arrange(desc(y)) %>%
    distinct(taxon, .keep_all = TRUE) %>%
    transmute(taxon, x = label_x, y)

  boxes <- tips %>%
    filter(!is.na(.data[[tax_rank]])) %>%
    mutate(taxon = .data[[tax_rank]]) %>%
    group_by(taxon) %>%
    summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")

  grad_layers <- lapply(seq_len(nrow(boxes)), function(i) {
    tax <- boxes$taxon[i]
    tax_fill <- if (tax %in% names(color_map)) color_map[[tax]] else "grey80"
    tax_grad <- grid::linearGradient(
      colours = c(scales::alpha(tax_fill, 0), scales::alpha(tax_fill, 0.65)),
      x1 = 0, y1 = 0.5, x2 = 1, y2 = 0.5)
    annotation_custom(
      grob = grid::rectGrob(gp = grid::gpar(fill = tax_grad, col = NA)),
      xmin = x_root, xmax = label_x,
      ymin = boxes$ymin[i], ymax = boxes$ymax[i])
  })

  grad <- grad +
    grad_layers +
    geom_tippoint(aes(color = family), shape = 15, size = 1.5) +
    scale_color_manual(values = family_colors, guide = "none") +
    geom_text2(data = labels, mapping = aes(x = x, y = y, label = taxon),
               inherit.aes = FALSE, size = 2.6, hjust = 0, vjust = 0.5)

  # Map() so each layer captures its own column (a for-loop would not)
  offsets <- c(first_offset, rep(bar_offset, length(bar_cols) - 1))
  bar_layers <- Map(function(col, off) {
    geom_fruit(
      geom = geom_bar,
      mapping = aes(x = !!rlang::sym(col), y = label, fill = !!rlang::sym(col)),
      orientation = "y", stat = "identity",
      offset = off, pwidth = bar_pwidth, width = bar_width)
  }, bar_cols, offsets)

  grad <- grad +
    bar_layers +
    scale_fill_fermenter(palette = "RdBu", n.breaks = 7,
                         limits = c(-bar_lim, bar_lim), name = bar_legend) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "right",
          plot.margin = margin(5.5, 90, 5.5, 5.5, "pt"))

  is_box <- vapply(grad$layers,
                   function(l) inherits(l$geom, "GeomCustomAnn"), logical(1))
  grad$layers <- c(grad$layers[is_box], grad$layers[!is_box])
  grad
}

# build_gift_tree_heatmap_len()
# Gradient tree + GIFT element heatmap of the degradation functions (D01 lipid, D02 polysaccharide, D03 sugar) + genome length bar. 
# Each element is centred and scaled across genomes.
build_gift_tree_heatmap_len <- function(grad_tree, gifts_sel,
                                        genome_metadata = macromag_genomemetadata) {
  gift_funcs_of_interest <- c("D01", "D02", "D03")

  gift_elements_long <- gifts_sel %>%
    to.elements(., GIFT_db) %>%
    as.data.frame() %>%
    rownames_to_column(var = "genome") %>%
    pivot_longer(cols = -genome, names_to = "Code_element", values_to = "GIFT") %>%
    filter(substr(Code_element, 1, 3) %in% gift_funcs_of_interest)

  uniqueGIFT_db <- unique(GIFT_db[, c("Code_element", "Code_function",
                                      "Function", "Element")])

  func_labels <- uniqueGIFT_db %>%
    distinct(Code_function, Function) %>%
    filter(Code_function %in% gift_funcs_of_interest) %>%
    mutate(Code_function = factor(Code_function, levels = gift_funcs_of_interest)) %>%
    arrange(Code_function)

  tip_y <- grad_tree$data %>%
    filter(isTip) %>%
    transmute(genome = label, y)

  # columns by function, then element code, with readable element names on top
  gift_heat_cols <- gift_elements_long %>%
    distinct(Code_element) %>%
    mutate(Code_function = factor(substr(Code_element, 1, 3),
                                  levels = gift_funcs_of_interest)) %>%
    left_join(uniqueGIFT_db %>% distinct(Code_element, Element), by = "Code_element") %>%
    arrange(Code_function, Code_element)
  gift_heat_col_levels <- gift_heat_cols$Code_element
  gift_heat_elem_labels <- setNames(gift_heat_cols$Element, gift_heat_cols$Code_element)

  gift_heat_mat <- gift_elements_long %>%
    filter(genome %in% tip_y$genome) %>%
    select(genome, Code_element, GIFT) %>%
    pivot_wider(names_from = Code_element, values_from = GIFT) %>%
    column_to_rownames("genome") %>%
    as.matrix()
  gift_heat_mat <- gift_heat_mat[, gift_heat_col_levels, drop = FALSE]

  gift_heat_z <- scale(gift_heat_mat)
  gift_heat_z[is.nan(gift_heat_z)] <- 0 # elements with sd = 0 give NaN
  gift_heat_values <- gift_heat_z
  gift_heat_zlim <- c(-ceiling(max(abs(gift_heat_z), na.rm = TRUE)),
                      ceiling(max(abs(gift_heat_z), na.rm = TRUE)))
  gift_heat_fill <- scale_fill_gradientn(
    colours = c("#7F5F00", "#D9A521", "white", "#2F7D88", "#0B3D46"),
    values = scales::rescale(c(gift_heat_zlim[1], gift_heat_zlim[1] / 2, 0,
                               gift_heat_zlim[2] / 2, gift_heat_zlim[2])),
    limits = gift_heat_zlim,
    oob = scales::squish,
    name = "GIFT (z)",
    guide = guide_colorbar(title.position = "top",
                           barwidth = unit(5, "cm"), barheight = unit(0.35, "cm"),
                           ticks.colour = "grey30"))

  gift_heat_long <- gift_heat_values %>%
    as.data.frame() %>%
    rownames_to_column("genome") %>%
    pivot_longer(-genome, names_to = "Code_element", values_to = "value") %>%
    mutate(Code_function = factor(substr(Code_element, 1, 3),
                                  levels = gift_funcs_of_interest),
           Code_element = factor(Code_element, levels = gift_heat_col_levels))

  gift_func_name_vec <- setNames(
    str_wrap(as.character(func_labels$Function), 16),
    as.character(func_labels$Code_function))

  gift_func_strip_colors <- c(D01 = "#F4D58D", D02 = "#BFD7B5", D03 = "#A9C5E0")

  gift_heat_func_facet <- ggh4x::facet_grid2(
    . ~ Code_function, scales = "free_x", space = "free_x",
    labeller = labeller(Code_function = gift_func_name_vec),
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = unname(gift_func_strip_colors[gift_funcs_of_interest]))))

  gift_heat_ylim <- c(min(tip_y$y) - 0.5, max(tip_y$y) + 0.5)

  gift_heat_tree <- grad_tree +
    scale_y_continuous(limits = gift_heat_ylim, expand = c(0, 0)) +
    coord_cartesian(clip = "off", ylim = gift_heat_ylim) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title.position = "top",
          plot.margin = margin(5.5, 70, 5.5, 5.5, "pt"))

  heat_panel <- gift_heat_long %>%
    inner_join(tip_y, by = "genome") %>%
    ggplot(aes(x = Code_element, y = y, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    gift_heat_fill +
    gift_heat_func_facet +
    scale_x_discrete(position = "top", labels = gift_heat_elem_labels) +
    scale_y_continuous(limits = gift_heat_ylim, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0, "pt"),
          axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank(), panel.spacing.x = unit(4, "pt"),
          strip.text = element_text(size = 7, face = "bold"),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.title.position = "top")

  len_panel <- tip_y %>%
    left_join(genome_metadata %>% select(genome, length), by = "genome") %>%
    mutate(panel = "Genome length") %>%
    ggplot(aes(x = "len", y = y, fill = length)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradientn(
      colours = c("#eef3f6", "#c9d7e1", "#98b2c2", "#66879c", "#35556b"),
      name = "Genome length",
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      guide = guide_colorbar(title.position = "top",
                             barwidth = unit(3, "cm"), barheight = unit(0.35, "cm"),
                             ticks.colour = "grey30")) +
    ggh4x::facet_grid2(. ~ panel,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(fill = "#c9d7e1"))) +
    scale_x_discrete(position = "top", labels = NULL) +
    scale_y_continuous(limits = gift_heat_ylim, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 6, "pt"),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          strip.text = element_text(size = 7, face = "bold"),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.title.position = "top")

  (gift_heat_tree + heat_panel + len_panel) +
    plot_layout(widths = c(0.7, 1.1, 0.10), guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.box = "vertical", legend.title.position = "top")
}

# build_extRLQ_scatter()
# An RLQ sample-side score against an environmental variable, with a linear trend and the Spearman correlation in the subtitle.
build_extRLQ_scatter <- function(df, x_col, axis_col, x_label, axis_label,
                                 axis_lim = NULL, legend_label = NULL,
                                 show_legend = TRUE) {
  ct <- suppressWarnings(cor.test(df[[x_col]], df[[axis_col]], method = "spearman"))
  lim <- if (is.null(axis_lim)) max(abs(df[[axis_col]]), na.rm = TRUE) else axis_lim
  leg_name <- if (is.null(legend_label)) axis_label else legend_label
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[axis_col]],
                      color = .data[[axis_col]])) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.5) +
    extRLQ_rdbu_scale(lim, leg_name) +
    labs(x = x_label, y = axis_label,
         subtitle = sprintf("Spearman rho = %.2f, p = %.3g",
                            unname(ct$estimate), ct$p.value)) +
    theme_minimal()
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# layout_extRLQ_scatter_pair()
# Two scatters side by side, y-axis on the left panel and legend on the right
layout_extRLQ_scatter_pair <- function(plot_left, plot_right) {
  plot_left <- plot_left + theme(legend.position = "none")
  plot_right <- plot_right +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  plot_left + plot_right + plot_layout(widths = c(1, 1))
}

# layout_extRLQ_scatter_grid()
# One column per cryosection, `plots_div` on top and `plots_seq` below
layout_extRLQ_scatter_grid <- function(plots_div, plots_seq) {
  n <- length(plots_div)
  stopifnot(length(plots_seq) == n)

  style_panel <- function(p, col, show_legend = FALSE) {
    p <- p + theme(legend.position = if (show_legend) "right" else "none")
    if (col > 1) {
      p <- p + theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.x = element_blank())
    }
    p
  }

  top <- Reduce(`|`, Map(function(p, i) style_panel(p, i), plots_div, seq_len(n)))
  bot <- Reduce(`|`, Map(function(p, i) style_panel(p, i, show_legend = (i == n)),
                         plots_seq, seq_len(n)))
  top / bot + plot_layout(widths = rep(1, n), heights = c(1, 1))
}


## Cryosection overlays

# cryosection_image()
cryosection_image <- function(cryo, dir = pixel_coords_dir, brightness = 100) {
  if (!dir.exists(dir)) return(NULL)
  hits <- list.files(dir, pattern = paste0("^", cryo, "_bright[.]png$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) {
    hits <- list.files(dir, pattern = paste0("^", cryo, "[.]png$"),
                       recursive = TRUE, full.names = TRUE)
  }
  if (length(hits) == 0) return(NULL)
  magick::image_modulate(magick::image_read(hits[1]), brightness = brightness)
}

# make_overlay_df()
# Pixel position of each microsample plus the value to colour it by. Image pixel y starts at the top, ggplot's at the bottom, so y is flipped
make_overlay_df <- function(cryo, space_df, value_col,
                            pixel_coords = pixel_coords_all, image_size = 1000) {
  pixel_coords %>%
    filter(cryosection == cryo, !is.na(pixel_x), !is.na(pixel_y)) %>%
    select(microsample, pixel_x, pixel_y) %>%
    inner_join(space_df %>% select(microsample, all_of(value_col)),
               by = "microsample") %>%
    mutate(pixel_y_flip = image_size - pixel_y)
}

# build_pixel_overlay()
# Microsamples drawn on cryosection image
build_pixel_overlay <- function(overlay_df, value_col, img, value_label,
                                diverging = TRUE, n_bins = 5, image_size = 1000,
                                axis_lim = NULL, show_legend = TRUE) {
  plot_df <- overlay_df
  if (diverging) {
    color_aes <- aes(x = pixel_x, y = pixel_y_flip, color = .data[[value_col]])
    lim <- if (is.null(axis_lim)) {
      max(abs(plot_df[[value_col]]), na.rm = TRUE)
    } else {
      axis_lim
    }
    color_scale <- extRLQ_rdbu_scale(lim, value_label)
  } else {
    bin_col <- paste0(value_col, "_bin")
    breaks <- unique(quantile(plot_df[[value_col]],
                              probs = seq(0, 1, length.out = n_bins + 1),
                              na.rm = TRUE))
    if (length(breaks) < 3) {
      breaks <- pretty(range(plot_df[[value_col]], na.rm = TRUE), n = n_bins)
    }
    plot_df[[bin_col]] <- cut(plot_df[[value_col]], breaks = breaks,
                              include.lowest = TRUE, dig.lab = 4)
    color_aes <- aes(x = pixel_x, y = pixel_y_flip, color = .data[[bin_col]])
    color_scale <- scale_color_viridis_d(name = value_label, option = "viridis")
  }

  ggplot(plot_df, color_aes) +
    annotation_raster(as.raster(img), xmin = 0, xmax = image_size,
                      ymin = 0, ymax = image_size) +
    geom_point(size = 2) +
    color_scale +
    coord_fixed(xlim = c(0, image_size), ylim = c(0, image_size), expand = FALSE) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
          panel.grid = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), plot.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          legend.position = if (show_legend) "right" else "none")
}

# plot_pixel_overlay()
# Overlay of one cryosection
plot_pixel_overlay <- function(cryo, space_df, value_col, value_label,
                               diverging = TRUE, axis_lim = NULL,
                               shift_x = 0, shift_y = 0,
                               pixel_coords = pixel_coords_all,
                               dir = pixel_coords_dir, ...) {
  img <- cryosection_image(cryo, dir = dir)
  if (is.null(pixel_coords) || is.null(img)) {
    message("no image or pixel coordinates for ", cryo, ", overlay skipped")
    return(invisible(NULL))
  }
  overlay_df <- make_overlay_df(cryo, space_df, value_col,
                                pixel_coords = pixel_coords) %>%
    mutate(pixel_x = pixel_x + shift_x, pixel_y_flip = pixel_y_flip + shift_y)
  build_pixel_overlay(overlay_df, value_col, img, value_label,
                      diverging = diverging, axis_lim = axis_lim, ...)
}

# show_overlay()
# Draw an overlay, or nothing when it was skipped
show_overlay <- function(p) if (is.null(p)) invisible(NULL) else print(p)
