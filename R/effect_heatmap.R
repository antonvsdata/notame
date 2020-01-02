
# Combines two matrices
coalesce<-function(...) {
  x<-lapply(list(...), function(z) {if (is.factor(z)) as.character(z) else z})
  m<-is.na(x[[1]])
  i<-2
  while(any(m) & i<=length(x)) {
    if ( length(x[[i]])==length(x[[1]])) {
      x[[1]][m]<-x[[i]][m]
    } else if (length(x[[i]])==1) {
      x[[1]][m]<-x[[i]]
    } else {
      stop(paste("length mismatch in argument",i," - found:", length( x[[i]] ),"expected:",length( x[[1]] ) ))
    }
    m<-is.na(x[[1]])
    i<-i+1
  }
  return(x[[1]])
}


#' Draw a heatmap of effects between variables, such as correlations
#'
#' Draws a heatmap of e.g. correlations between variables (see perform_correlation_tests).
#' It is possible to draw only the lower triangular of the heatmap, order rows and columns
#' with hierarchical clustering, and add circles for p-values.
#'
#' @param data a data frame with x and y variables and the effect (used to fill the tiles)
#' @param x,y the column names of data with the x and y variables
#' @param effect the column name of the effect, e.g. correlation
#' @param p optional, the column name with p-values. If provided, points that scale by p-value are drawn
#' on top of the heatmap tiles
#' @param p_limit numeric, only p-values below the limit are plotted as points
#' @param point_size_range a numeric vector of length 2. The upper and lower limits for the point sizes.
#' This needs to be adjusted to make the point size look good when compared to the tiles
#' @param log2_effect logical, whether the effect should be plotted on a logarithmic scale (in case of fold change etc.)
#' @param discretize_effect logical, whether the effect range should be divided into discrete levels instead of using
#' a continuous scale. Can sometimes make patterns more visible, but the hard limits can blur the big picture as well.
#' @param breaks if \code{discretize_effect = TRUE}, either the number of breaks or the points where to cut for the levels,
#' see \code{\link{cut}}
#' @param clustering logical, whether the order of rows and columns should be ordered by hierarchical clustering?
#' @param dist_method distance method used in clustering, see \code{\link{dist}}
#' @param clust_method clustering method used in clustering, see \code{\link{hclust}}
#' @param lower_tri logical, should only the lower triangular be plotted?
#' @param reverse_y logical, if \code{clustering = FALSE, lower_tri = FALSE}, should the order of the y-axis
#' be reversed so that the diagonal is from top left to bottom right?
#' @param title,subtitle the title and subtitle of the plot
#' @param fill_scale fill scale for the heatmap as returned by a ggplot function. Set to NA to choose the appropriate scale based on the class of the effect variable.
#'
#' @return a ggplot object
#'
#' @details All missing effects between variables are replaced by 0 before clustering,
#' since \code{hclust} can't deal with missing values.
#'
#' @seealso \code{\link{cut}} for discretizing the effect, \code{\link{dist}} for distance calculation for clustering,
#' \code{\link{hclust}} for hierarchical clustering.
#'
#' @examples
#' # Compute correlations between variables
#' correlations <- perform_correlation_tests(example_set, x = featureNames(example_set), duplicates = TRUE)
#'
#' # Minimal example
#' plot_effect_heatmap(correlations, x = "X", y = "Y", effect = "Correlation_coefficient")
#'
#' # Lower triangular with discrete effect and p-value dots
#' plot_effect_heatmap(correlations, x = "X", y = "Y", effect = "Correlation_coefficient", p = "Correlation_P",
#'                point_size_range = c(2,8),
#'                discretize_effect = TRUE, breaks = 7, lower_tri = TRUE)
#'
#' @export
plot_effect_heatmap <- function(data, x, y, effect, p = NULL, p_limit = 0.1, point_size_range = c(1, 6),
                           log2_effect = FALSE, discretize_effect = FALSE, breaks = 5,
                           clustering = TRUE, dist_method = "euclidean", clust_method = "ward.D2",
                           lower_tri = FALSE, reverse_y = TRUE,
                           title = NULL, subtitle = NULL, fill_scale = NA) {
  # Get default fill scales
  if (!is.null(fill_scale)) {
    if (is.na(fill_scale)) {
      if (discretize_effect | class(data[, effect]) %in% c("factor", "character")) {
        fill_scale <- getOption("notame.fill_scale_div_dis")
      } else {
        fill_scale <- getOption("notame.fill_scale_div_con")
      }
    }
  }

  # Possible log-transform effect, should show on legend
  legend_label <- effect
  if (log2_effect) {
    data[, effect] <- log2(data[, effect])
    legend_label <- paste0("log2(", effect, ")")
  }

  if (lower_tri) {
    # Clustering is handled inside to_lowertri
    data <- to_lowertri(data, x, y, effect, clustering, clust_method, dist_method)
  } else if (clustering) {
    data <- hclust_effects(data, x, y, effect, clust_method, dist_method)
  } else if(reverse_y) {
    # Reverse order of y axis so diagonal is from top left to bottom right
    data[y] <- factor(data[, y], levels = rev(levels(factor(data[, y]))))
  }

  if (discretize_effect) {
    data[effect] <- cut(data[, effect], breaks = breaks,
                        dig.lab = 2)
    data[effect] <- factor(data[, effect], levels = rev(levels(data[, effect])))
  }

  ggp <- ggplot(data, aes_string(x = x , y = y, fill = effect)) +
    geom_tile() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    labs(x = "", y = "", fill = legend_label,
         title = title, subtitle = subtitle) +
    coord_fixed() +
    theme(aspect.ratio=1) +
    fill_scale

  if (!is.null(p)) {
    small_p <- data[data[,p] < p_limit, ]
    small_p[p] <- -log10(small_p[, p])
    ggp <- ggp +
      geom_point(aes_string(size = p), data = small_p,
                 colour = 'grey10', fill = 'grey30',alpha = 0.3, shape = 21) +
      scale_size(name = '-log10(p-value)', range = point_size_range)
  }

  ggp
}

# Helper function to cluster effects without setting to lower triangular
hclust_effects <- function(data, x, y, effect, clust_method, dist_method) {

  # Convert to wide format matrix for clustering
  data_wide <- data %>%
    dplyr::select(x, y, effect) %>%
    tidyr::spread(y, effect) %>%
    dplyr::filter(!is.na(x)) %>%
    tibble::column_to_rownames(x) %>%
    as.matrix()

  data_wide[is.na(data_wide)] <- 0

  # Separate clustering for y and x
  x_order <- hclust(dist(data_wide, method = dist_method), method = clust_method)$order
  y_order <- hclust(dist(t(data_wide), method = dist_method), method = clust_method)$order
  data_wide <- data_wide[x_order, y_order]

  data[x] <- factor(data[, x], levels = rev(rownames(data_wide)))
  data[y] <- factor(data[, y], levels = colnames(data_wide))
  data
}



# Converts data to only include the lower triangular
to_lowertri <- function(data, x, y, effect, clustering, clust_method, dist_method) {


  # Rename columns for simplicity
  data <- data %>% dplyr::rename("x" = x, "y" = y, "effect" = effect)

  dat <- data[c("x", "y", "effect")]
  dat$x <- as.character(dat$x)
  dat$y <- as.character(dat$y)

  # This makes sure all the variables will be included in both axes
  vars <- unique(c(dat$x, dat$y))
  x_missing <- setdiff(vars, dat$x)
  y_missing <- setdiff(vars, dat$y)
  append_len <- max(length(x_missing),length(y_missing))
  if (append_len) {
    append_df <- data.frame(x = c(x_missing,rep(NA, append_len - length(x_missing))),
                            y = c(y_missing,rep(NA, append_len - length(y_missing))),
                            effect = rep(NA, append_len))
    colnames(append_df) <- colnames(dat)

    dat <- rbind(dat,append_df)
  }


  # Converting data into wide format and tidying data
  dat_w <- dat %>% tidyr::spread(y, effect) %>% dplyr::filter(!is.na(x))
  # Move column x to rownames
  dat_w <- tibble::column_to_rownames(dat_w, "x")
  dat_w <- dat_w[rownames(dat_w)!="NA",! colnames(dat_w) %in% c("NA","<NA>")] #one of the columns or rows is named NA
  dat_w <- dat_w[rev(order(names(dat_w))),rev(order(names(dat_w)))] #this makes the image lie on the lower triangular


  # dat_w only has one-directional interactions
  # for example, metformine vs morphine = 1.09, but morphine vs metformine = NA
  # transpose values are added to make the matrix symmetrical
  dat_w_t <- data.frame(t(dat_w))
  dat_w_whole <- coalesce(dat_w,dat_w_t)

  if (clustering){
    data_w_zeros <- dat_w_whole
    data_w_zeros[is.na(data_w_zeros)] <- 0
    hc <- hclust(dist(data_w_zeros, dist_method), method = clust_method)
    dat_w_whole <- dat_w_whole[hc$order, hc$order]
    dat_w_whole
  }

  # Only half of the associations are needed for plotting
  dat_w_whole[upper.tri(dat_w_whole)] <- NA
  # Diagonal should be included in the plot


  # Melt back to long format for ggplot2
  dat_w_whole$x <- rownames(dat_w_whole)
  dat_l <- tidyr::gather(dat_w_whole, y, effect, -x)
  dat_l$y <- as.character(dat_l$y)


  # The order of Variable1 and 2 has changed for some associations, so two joins are required
  combined1 <- dplyr::inner_join(dat_l, data, by = c("x", "y", "effect"))
  combined2 <- dplyr::inner_join(dat_l, data, by = c("x" = "y","y" = "x","effect"))
  dat_l <- rbind(combined1,combined2)%>%
    dplyr::distinct() # Remove duplicated associations with same x and y

  # Setting the factor levels to correctly draw the heatmap
  # This ensures the tiles are plotted in correct order to make a lower triangular heat map
  dat_l$x <- dat_l$x %>%
    factor(levels = rev(rownames(dat_w_whole)))
  dat_l$y <- dat_l$y %>%
    factor(levels = rownames(dat_w_whole))

  colnames(dat_l)[colnames(dat_l) == "x"] <- x
  colnames(dat_l)[colnames(dat_l) == "y"] <- y
  colnames(dat_l)[colnames(dat_l) == "effect"] <- effect

  dat_l

}
