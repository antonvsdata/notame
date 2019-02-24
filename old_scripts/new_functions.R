library(memoise)

# The srawing function used by save.ordered.boxplot
draw.ordered.boxplot <- function(data, responses, group.col, order.col,
                              zoom.boxplot=TRUE, title="Sample boxplot from metabolite values") {

  require(reshape)
  # Generate new sample id
  data$QC <- data[,group.col] == "QC"
  #data <- na.omit(data[c("QC", group.col, order.col, responses)])

  # Melt data for ggplot to make boxplots from multiple responses
  data <- melt(data,
               id.vars=c("QC", group.col, order.col),
               measure.vars=responses)

  # Order data by order.col
  data <- data[order(data[,order.col]),]

  # Make sure x.levels as factor and get levels for ordering
  x.levels <- levels(factor(data[,order.col], levels=unique(as.character(data[,order.col]))))

  data[,order.col] <- as.factor(as.character(data[,order.col]))

  # Plotting

  p <- ggplot(data, aes_string(x=order.col, y="value"))#, group=group.col, fill=group.col))
  # Colour boxplots by Group
  p <- p + geom_boxplot(aes_string(fill="QC")) +
    scale_fill_manual(values=c("#33ccff","#007399"))

  ## Zooming outliers out of view
  if(zoom.boxplot){
    # compute lower and upper whiskers
    ylim1 = boxplot.stats(data$value)$stats[c(1, 5)]
    # scale y limits based on ylim1
    p <- p + coord_cartesian(ylim = ylim1*1.05)
    # add text to main title
    title <- paste(title, "(zoomed in boxplot: outliers out of view)", sep="")
  }

  # Get specific x-axis order for the plot
  p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)

  # Add title, apply theme
  p <- p + theme(legend.position="top")
  p <- p + ggtitle(title) +
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, hjust=0.5))
  p <- p + theme_bw()

  p <- p + theme(axis.text.x = element_text(angle=90, hjust=0.5))

  p
}

# Save boxplots of all response values, one boxplot per sample, ordered by values in order.col, e.g. injection order
# group.col should ba a character or factor, and have the value "QC" for QC samples and other values for biological samples
# The boxplots are colored with a darker shade for QC samples
# Set zoom.boxplot = TRUE to zoom in and exclude outliers
save.ordered.box.plot <- function(data, responses, group.col, order.col, zoom.boxplot=FALSE, output.file, log.file, logging){
  # Create plot
  p <- draw.ordered.boxplot(data=data, responses=responses,
                         group.col = group.col,
                         order.col = order.col,
                         zoom.boxplot=TRUE)
  # ggSave plot, logging in ggsave.plot
  ggsave.plot(p, plot.file=output.file, log.description="Saved ordered samples boxplots - plot to: ", log.file, logging, scale=2.5)
}

# PCA plot with continuous coloring by one variable
# color.col is the coloring column name as a character
# color.scale:  "blues"   continuous blue color scale
#               "grey"    continuous greyscale
#               "3color"  blue-white-red color scale, specify midpoint, defaults to mean
save.pca.colored <- function(data, responses, color.col, color.scale = "blues", color.midpoint = NULL, width=8, height=6, output.file, log.file, logging){
  require(gridExtra)

  pca.result <- compute.pca(data[responses])

  # Same PCA plot without sample IDs
  p <- plot.pca.colored(pca.data=pca.result, color.ref = data[,color.col], color.title= color.col, color.scale = color.scale, color.midpoint = color.midpoint)

  #ggsave.plot(p, plot.file=output.file, log.description="Saved PCA plot to: ", log.file, logging, scale=1.5, width=6, height=12)
  save.plot(p, plot.file=output.file, log.description="Saved colored PCA plot to: ", log.file=log.file, logging=logging, width=width, height=height)
}

# Drawing function used by save.pca.colored
plot.pca.colored <- function(pca.data, color.ref, color.title, color.scale, color.midpoint, component.x = 1, component.y = 2){

  require(MASS)
  require(ggplot2)
  require(scales)
  #require(gridExtra)

  # Capture local enviroment, because problematic ggplot aes with local variables
  .e <- environment()

  # Grouping variable and pca results to same dataframe
  scores <- data.frame(coloring=color.ref, pca=pca.data$x)
  scores$coloring <- as.numeric(as.character(scores$coloring))


  # for percentage labels
  prop.pca = pca.data$sdev^2/sum(pca.data$sdev^2)

  # string variables made from component numbers to be referenced
  PCX <- paste("pca.PC", as.character(component.x), sep="")
  PCY <- paste("pca.PC", as.character(component.y), sep="")


  p <- ggplot(scores, environment=.e) +
    geom_point(aes_string(x=PCX, y=PCY, colour="coloring"),size=2.5) +
    labs(x = paste(PCX, "(", percent(prop.pca[component.x]), ")", sep=""),
         y = paste(PCY, "(", percent(prop.pca[component.y]), ")", sep=""),
         color=color.title,
         title = "PCA") +
    theme_bw()

  if(color.scale == "grey" | color.scale == "gray"){
    p <- p +
      scale_colour_gradient(low = "grey80", high = "grey20")
  }
  if(color.scale == "3colors"){
    if(is.null(color.midpoint)){
      color.midpoint <- mean(scores$coloring)
    }
    p <- p +
      scale_colour_gradient2(low = "steelblue", mid = "white", high = "red", midpoint = color.midpoint, space = "Lab", na.value = "green")
  }
  p
}

# Save tSNE plot colored by injection order (or another coloring variable)
# oreder.col is the ordering column name as character
save.tsne.plot.order <- function(data, responses, order.col, output.file, logging, log.file){
  require(Rtsne)

  tmp.matrix <- as.matrix(data[responses])
  set.seed(42)
  data.rtsne <- Rtsne(tmp.matrix)

  #plot(dat.rtsne$Y, col=data$GROUP)

  # Dataframe for ggplot
  temp.df <- data.frame(data.rtsne$Y, ORDER=data[,order.col])
  # Change name of the (grouping) column
  colnames(temp.df) <- c("X1","X2",order.col)

  # Plot with ggplot2
  p <- ggplot(temp.df, aes_string(x="X1", y="X2", colour=order.col)) +
    geom_point()
  # add theme and title
  p <- p + theme(legend.position="top")
  p <- p + ggtitle("TSNE") +
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, hjust=0.5))
  p <- p + theme_bw()

  # X,Y-axis to same aspect ratio
  p <- p + coord_fixed()

  # Save plot
  pdf(output.file)
  plot(p)
  dev.off()

  log.text(paste("Saved tSNE plot to: ", output.file, sep=""), log.file=log.file, logging=logging)
}

# t-SNE plot with subjecct IDs on another page
save.tsne.plot.id <- function(data, responses, group.col, id.col = NULL, output.file, width = 10, height = 10, logging, log.file){
  require(Rtsne)
  if(!is.null(id.col)){
    require(ggrepel)
  }
  
  tmp.matrix <- as.matrix(data[responses])
  set.seed(42)
  data.rtsne <- Rtsne(tmp.matrix)
  
  #plot(dat.rtsne$Y, col=data$GROUP)
  
  # Dataframe for ggplot
  if(!is.null(id.col)){
    temp.df <- data.frame(data.rtsne$Y, GROUP=data[,group.col], ID = data[,id.col])
  } else {
    temp.df <- data.frame(data.rtsne$Y, GROUP=data[,group.col])
  }
  
  # Change name of the (grouping) column
  colnames(temp.df) <- c("X1","X2", group.col, "ID")
  
  # Plot with ggplot2
  p <- ggplot(temp.df, aes_string(x="X1", y="X2")) +
    geom_point(aes_string(colour=group.col, shape=group.col, group=group.col))
  # add theme and title
  p <- p + theme(legend.position="top")
  p <- p + ggtitle("TSNE") + 
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, hjust=0.5))
  p <- p + theme_bw()
  
  # X,Y-axis to same aspect ratio
  p <- p + coord_fixed()
  
  # Save plot
  pdf(output.file, width = width, height = height, onefile = TRUE)
  plot(p)
  if(!is.null(id.col)){
    p2 <- p +
      geom_text_repel(aes(label = ID))
    plot(p2)
  }
  dev.off()
  
  log.text(paste("Saved tSNE plot to: ", output.file, sep=""), log.file=log.file, logging=logging)
}

# Dendrogram plot of hierarchical clustering
save.dendrograms <- function(data, responses, id.col, group.col, apply.scale = TRUE, distance.metrics = "euclidean", clustering.methods = c("complete", "ward.D", "ward.D2"),
                             output.file, width = 50, height = 16, log.file, logging){
  require(dplyr)
  require(ggplot2)
  require(ggdendro)
  
  rownames(data) = data[, id.col]
  
  settings <- expand.grid(clust = clustering.methods, dist = distance.metrics, stringsAsFactors = FALSE)
  
  clust.data <- data[responses]
  if(apply.scale){
    clust.data <- scale(clust.data)
  }
  pdf(output.file, onefile = TRUE, width = width, height = height)
  for(i in 1:nrow(settings)){
    ddata <- hclust(dist(clust.data, method = settings$dist[i]), method = settings$clust[i]) %>%
      as.dendrogram() %>%
      dendro_data()
    
    labels <- label(ddata) %>%
      dplyr::mutate(label = as.character(label)) %>%
      dplyr::left_join(data[c(id.col, group.col)], by = c("label" = id.col))
    labels[, group.col] <- as.factor(labels[, group.col])
    
    p <- ggplot(segment(ddata)) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data = labels, aes_string(x="x", y="y", label = "label", color = group.col), angle = 90, hjust = 1) +
      theme_dendro() +
      ggtitle(paste("Ditance metric:", settings$dist[i], "  Clustering method:", settings$clust[i]))
    plot(p)
  }
  
  dev.off()
  log.text(paste("Saved dendrogram to:", output.file), log.file = log.file, logging = logging)
}


# Hexbin plot of PCA results, hexagons colored by one variable
# If the coloring variable is not given, hexagons are colored by count of samples in them
# color.col is the coloring column name as a character
# summary.function defines the summary function used when coloring column is specified
# color.scale:  "blues"   continuous blue color scale
#               "grey"    continuous greyscale
#               "3color"  blue-white-red color scale, specify midpoint, defaults to mean
save.pca.hexbin <- function(data, responses, color.col, summary.function = "mean", color.scale = "blues", color.midpoint = NULL, width=8, height=6, output.file, log.file, logging){
  require(gridExtra)

  pca.result <- compute.pca(data[responses])

  # Same PCA plot without sample IDs
  pca.plot <- plot.pca.hexbin(pca.data=pca.result, color.ref=data[,color.col], color.title=color.col, summary.function = summary.function, color.scale = color.scale, color.midpoint = color.midpoint)

  #ggsave.plot(p, plot.file=output.file, log.description="Saved PCA plot to: ", log.file, logging, scale=1.5, width=6, height=12)
  save.plot(pca.plot, plot.file=output.file, log.description="Saved hexbin PCA plot to: ", log.file=log.file, logging=logging, width=width, height=height)
}

# New function for plotting the hexbin plot of PCA
# Color scale: single or 3color
# For three colors, specify the midpoint
plot.pca.hexbin <- function(pca.data, color.ref, color.title, summary.function, color.scale, color.midpoint, component.x = 1, component.y = 2){

  require(MASS)
  require(ggplot2)
  require(scales)
  #require(gridExtra)

  # Capture local enviroment, because problematic ggplot aes with local variables
  .e <- environment()

  # Grouping variable and pca results to same dataframe

  scores <- data.frame(coloring=color.ref, pca=pca.data$x)
  scores$coloring <- as.numeric(as.character(scores$coloring))

  # for percentage labels
  prop.pca = pca.data$sdev^2/sum(pca.data$sdev^2)

  # string variables made from component numbers to be referenced
  PCX <- paste("pca.PC", as.character(component.x), sep="")
  PCY <- paste("pca.PC", as.character(component.y), sep="")
  p <- ggplot(scores, environment = .e) +
    stat_summary_hex(aes_string(x=PCX, y=PCY, z = "coloring"), bins = 10, fun = summary.function)+
    labs(x = paste(PCX, "(", percent(prop.pca[component.x]), ")", sep=""),
         y = paste(PCY, "(", percent(prop.pca[component.y]), ")", sep=""),
         fill=paste(summary.function,"(", color.title,")",sep=""),
         title="PCA")+
    theme_bw()
  if(color.scale == "grey"){
    p <- p +
      scale_fill_gradient(low = "grey80", high = "grey20")
  }
  if(color.scale == "3colors"){
    if(is.null(color.midpoint)){
      color.midpoint <- mean(scores$coloring)
    }
    p <- p +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "red", midpoint = color.midpoint, space = "Lab", na.value = "green")
  }

  p
}

# Compoute linear regression model response ~ phenotype for each response
perform.lm <- function(data, responses, phenotype, covariate.string=NULL) {

  results <- NULL
  results <- foreach(i=1:length(responses), .combine=rbind) %dopar% {
    response <- responses[i]
    if (is.null(covariate.string)) {
      formula <- as.formula(paste (response, " ~ ", phenotype, sep=""))
    } else {
      formula <- as.formula(paste (response, " ~ ", phenotype, covariate.string, sep=""))

    }
    fit <- NULL
    tryCatch({
      fit <- lm(formula, data=data)
    }, error = function(e) print(e$message))
    if(is.null(fit) | sum(!is.na(data[,response])) < 2){
      result.row <- data.frame(Trait=response, Rsq=NA, SD=NA, L95CI=NA, U95CI=NA, Estimate=NA, P=NA)
    } else {
      #print(summary(fit))
      phenotype_ind <- paste(phenotype,levels(data[,phenotype])[2],sep="")
      p <- summary(fit)$coefficients[phenotype_ind,"Pr(>|t|)"]
      sd <- summary(fit)$coefficients[phenotype_ind,"Std. Error"]
      est <- summary(fit)$coefficients[phenotype_ind,"Estimate"]
      rsq <- summary(fit)$r.squared
      l95 <- confint(fit, level=0.95)[phenotype_ind,1]
      u95 <- confint(fit, level=0.95)[phenotype_ind,2]

      result.row <- data.frame(Trait=response, Rsq=rsq, SD=sd, L95CI=l95, U95CI=u95, Estimate=est, P=p)
    }
    result.row
  }
  results$"P-FDR" <- p.adjust(results$P, method="BH")
  results
}

# Save histograms of p-values of linear regression coefficients between metabolite levels and values in lm.col
# group.col should contain groups with QCs, calculation will be run separately for QC samples
save.lm.pvalue.histograms <- function(data, responses, group.col, lm.col, bin.width = 0.05, output.file, log.file, logging){

  lm.results <- perform.lm(data=data, responses = responses, phenotype = lm.col)

  data.QCs <- data[data[group.col] == "QC",]

  lm.results.qc <- perform.lm(data=data.QCs, responses = responses, phenotype = lm.col)

  p.values <- data.frame(ALL_SAMPLES = lm.results$P, QC_SAMPLES = lm.results.qc$P)

  save.pvalue.histograms(p.values, output.file=output.file, log.file=log.file, logging=logging)
}

# used by save.lm.pvalue.histograms
save.pvalue.histograms <- function(p.values, bin.width = 0.05, width = 6, height = 6, output.file, log.file, logging){

  require(gridExtra)

  plots <- draw.pvalue.histograms(p.values, bin.width)

  p <- arrangeGrob(plots[[1]], plots[[2]], nrow = 2)

  save.plot(p, plot.file=output.file, log.description="Saved p-value histograms to: ", log.file=log.file, logging=logging, width=width, height=height)
}

# used by save.lm.pvalue.histograms
draw.pvalue.histograms <- function(p.values, bin.width){
  breaks <- seq(0,1,by=bin.width)


  plots <- foreach(i = 1:ncol(p.values), .combine = c, .packages = "ggplot2") %dopar% {
    p.finite <- p.values[,i][which(is.finite(p.values[,i]))]
    y.line <- length(p.finite)/(length(breaks)-1)

    colname <- colnames(p.values)[i]
    p <- ggplot(p.values) +
      geom_histogram(aes_string(colname), breaks = breaks, col = "grey20", fill = "white") +
      geom_hline(yintercept = y.line, color="red", linetype = "dashed") +
      labs(x="P-value", y="Frequency") +
      ggtitle(colname) +
      theme_bw() +
      theme(plot.title = element_text(face="bold", hjust=0.5))

    list(p)
  }

  plots
}

# Save density plots of distances between samples
# group.col should contain groups with QCs, calculation will be run separately for QC samples and non-QC samples
# specify the distance method used with dist.method, defaults to euclidean.
save.distance.density.plot <- function(data, responses, dist.method = "euclidean", group.col, width=6, height=6, output.file, log.file, logging){

  title <- paste("Density plot of",dist.method,"distances between samples")

  data.noQCs <- data[data[,group.col] != "QC",]
  D.noqc <- dist(data.noQCs[,responses], method = dist.method) %>% as.numeric()
  data.QCs <- data[data[,group.col] == "QC",]
  D.qc <- dist(data.QCs[,responses], method = dist.method) %>% as.numeric()
  Group <- c(rep("non-QC",length(D.noqc)), rep("QC",length(D.qc)))
  distances <- data.frame(dist = c(D.noqc,D.qc), Group = Group)

  p <- ggplot(distances) +
    geom_density(aes(dist, fill = Group, colour = Group), alpha = 0.2) +
    ggtitle(title) +
    theme_minimal()

  save.plot(p, plot.file=output.file, log.description="Saved density plot to: ", log.file=log.file, logging=logging, width=width, height=height)
}

# Choose number of principal components to use in the model for residuals plot
# If the percentage of variance explained by a component is lower than value.lmit AND
# the difference in percentages of variance explained to the next component is below diff.limit
# the component will be the last component included in the model
choose.n.comp <- function(eig, value.limit, diff.limit){
  # Convert eigenvalues to percentages of variance explained
  perc <- eig/sum(eig)
  dif <- abs(diff(perc))
  n.comp <- min(which(perc[-length(perc)] < value.limit & dif < diff.limit))
  if(n.comp == Inf){
    n.comp = length(eig)
  }
  n.comp
}

# Draw a Q-residuals vs Hotelling T2 plot of PCA
# n.comp (optional) specifies the number of components. if left NULL, the optimal number will be chosen according to
#   value.limit and diff.limit (see choose.n.comp above)
# group.col (optional) specifies column to use for coloring
# id.col (optional) specifies column to use for point labels
save.pca.residuals.plot <- function(data, responses, n.comp = NULL, value.limit = 0.03, diff.limit = 0.004, group.col = NULL, id.col = NULL, width=6, height=6, output.file, log.file, logging){

  require(mdatools)

  if(!is.null(group.col)){
    c.group <- data[,group.col]
    if(class(c.group) != "factor"){
      c.group <- as.factor(c.group)
    }
  }
  else{
    c.group <- NULL
  }
  if(!is.null(id.col) & length(unique(data[,id.col])) != nrow(data)){
    stop("All values in the ID column are not unique")
  }

  m <- mdatools::pca(data[responses])
  if(is.null(n.comp)){
    n.comp <- choose.n.comp(m$eigenvals, value.limit, diff.limit)
  }
  m <- selectCompNum(m,n.comp)

  output.file <- gsub(".pdf","",output.file)

  pdf(file=paste(output.file, n.comp,"PCs.pdf",sep="_"), width=width, height=height, onefile = TRUE)
  plotResiduals(m$calres, cgroup = c.group, show.labels = FALSE)
  # If id.col is specified, plot another PCA with sample IDs
  if(!is.null(id.col)){
    rownames(data) <- as.character(data[,id.col])
    m <- mdatools::pca(data[responses])
    m <- selectCompNum(m,n.comp)
    plotResiduals(m$calres, cgroup = c.group, show.labels = TRUE)
  }
  dev.off()
  log.text(paste("Saved PCA residuals plot with", n.comp, "principal components to:",output.file), log.file = log.file, logging = logging)
}

# Save plot of variance explained by first 15 principal components
save.pca.variance.plot <- function(data, responses, width=6, height=6, output.file, log.file, logging){

  require(mdatools)

  m <- mdatools::pca(data[responses])

  pdf(file=output.file, width=width, height=height)
  plotVariance(m, type='b')
  dev.off()
  log.text(paste("Saved PCA variance plot to:",output.file), log.file = log.file, logging = logging)
}

# Save a plot of PCA loadings
# response.labels = TRUE (default) shows response names in the plot
# point and text size as well as width and height can be set, the defaults should work OK
save.pca.loadings.plot <- function(data, responses, n.responses, output.file, plot.title = "PCA loadings", point.size = 2, text.size = 3, width = 10, height = 10, fixed_coordinates=FALSE,
                                   response.labels = TRUE, log.file, logging){
  pca.data <- compute.pca(data[responses])

  p <- plot.pca.loadings(pca.data = pca.data, n.responses = n.responses, title = plot.title, point.size = point.size,
                         text.size = text.size, fixed_coordinates = fixed_coordinates, response.labels = response.labels)

  pdf(output.file, width = width, height = height)
  print(p)
  dev.off()
  log.text(paste("Saved PCA loadings plot to:", output.file), log.file = log.file, logging = logging)
}

# used by save.pca.loadings plot
plot.pca.loadings <- function(pca.data, n.responses, title, point.size, text.size, fixed_coordinates, response.labels){

  require(ggrepel)
  require(scales)
  pca <- data.frame(response = rownames(pca.data$rotation), pca.data$rotation[,1:2],
                    importance = abs(pca.data$rotation[,1] * pca.data$sdev[1]^2) + abs(pca.data$rotation[,2] * pca.data$sdev[2]^2))
  pca <- pca[order(pca$importance, decreasing = TRUE)[1:n.responses],]

  prop.pca <- pca.data$sdev^2/sum(pca.data$sdev^2)

  p <- ggplot(pca, aes(x = PC1, y = PC2)) +
    geom_point(size = point.size) +
    labs(x = paste0("PC1 (", percent(prop.pca[1]), ")"),
         y = paste0("PC2 (", percent(prop.pca[2]), ")")) +
    theme_minimal() +
    ggtitle(title) +
    geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey80", linetype = "dashed")

  if(fixed_coordinates){
    p <- p + coord_fixed()
  }
  if(response.labels){
    p <- p +
      geom_text_repel(aes(label = response), size = text.size)
  }

  p

}

# Save histograms of different quality measurements
# data = result from compute quality measurements
# cutoffs is a numeric vector of length 4, giving the cutoff values for each quality measure
# the order is: RSD, RSD_r, D_ratio, D_ratio_r
save.column.histograms <- function(data, plot.file, columns = NULL, cutoffs = NULL, bin.widths = NULL, log.file, logging){
  # If columns are not given, default to all numeric columns
  if(is.null(columns)){
    col.classes <- sapply(data, class)
    columns <- colnames(data)[col.classes == "numeric"]
  }

  plots <- foreach(i = 1:length(columns), .combine = c, .packages = "ggplot2") %do% {
    colname = columns[i]
    filt <- (is.finite(data[,colname])) & (abs(data[,colname]) < 10*IQR(data[,colname], na.rm = TRUE))
    data.tmp <- data[filt,]
    values <- data.tmp[,colname]
    # If all the values are the same (like detection rate in impted data), skip the column
    if(is.na(var(values, na.rm = TRUE)) | var(values, na.rm = TRUE) == 0){
      return(NULL)
    }
    # If the bin.width is not given, compute it with the Freedman-Diaconis rule
    if(is.null(bin.widths)){
      bin.width <- 2 * IQR(values) / length(values)^(1/3)
    }else{
      bin.width <- bin.widths[i]
    }
    # Create the histogram
    breaks <- seq(min(values), max(values), by = bin.width)
    p <- ggplot(data.tmp) +
      geom_histogram(aes_string(colname), breaks = breaks, col = "grey20", fill = "white") +
      labs(x=colname, y="Frequency") +
      ggtitle(colname) +
      theme_bw() +
      theme(plot.title = element_text(face="bold", hjust=0.5))
    # If cutoffs are specified, draw red vertical line at cutoff value
    if(!is.null(cutoffs)){
      p <- p + geom_vline(xintercept = cutoffs[i], color="red", linetype = "dashed")
    }
    list(p)
  }
  # Draw all the plots as separated pages in the same pdf
  pdf(plot.file, onefile=TRUE)
  for (x in 1:length(plots)){
    print(plots[[x]])
  }
  dev.off()
  log.text(paste("Save quality measurement histograms to:", plot.file), log.file = log.file, logging = logging)
}

draw.dot.errorbar.plot <- function(data, response, group.col, subtitle, errorbar.fun,
             errorbar.centre, errorbar.min, errorbar.max, main.plot = FALSE){
  p <- ggplot(data, aes_string(x=group.col, y=response, colour=group.col)) +
    # Errorbars with solid lines
    stat_summary(fun.data =errorbar.fun,
                 geom="errorbar", width = 1,
                 fun.y = errorbar.centre, fun.ymin = errorbar.min, fun.ymax = errorbar.max) +
    # Plot point to mean
    stat_summary(fun.data = mean_se,
                 geom = "point",
                 size = 4) +
    theme_bw() +
    scale_colour_brewer(type = "qual", palette = "Dark2")
  if(main.plot){
    p <- p +
      labs(title = response, subtitle = subtitle)
  } else {
    p <- p +
      labs(subtitle = subtitle)
  }

  p
}

save.dot.errorbar.plots <- function(data1, data2 = NULL, responses, group.col, data1.title = NULL, data2.title = NULL, errorbar.fun = "mean_cl_boot",
                          errorbar.centre = NULL, errorbar.min = NULL, errorbar.max = NULL, output.file, log.file, logging){
  require(ggplot2)
  require(gridExtra)
  require(Hmisc)
  
  # Strip NAs
  data1 <- data1[!is.na(data1[,group.col]), ]
  data2 <- data1[!is.na(data2[,group.col]), ]
  
  pdf(output.file, onefile = TRUE)
  for(response in responses){

    p <- draw.dot.errorbar.plot(data1, response, group.col, subtitle = data1.title, errorbar.fun,
                      errorbar.centre, errorbar.min, errorbar.max, main.plot = TRUE)

    if(!is.null(data2)){
      p2 <- draw.dot.errorbar.plot(data2, response, group.col, subtitle = data2.title, errorbar.fun,
                         errorbar.centre, errorbar.min, errorbar.max)
      p <- arrangeGrob(p, p2, nrow = 2)
    }
    plot(p)
  }
  dev.off()
  log.text(paste("Saved dot plots with errorbars to:", output.file), log.file = log.file, logging = logging)
}

# Helper functions for quality measurement computations
finite.sd <- function(x){
  sd(x[is.finite(x)], na.rm = TRUE)
}

finite.mean <- function(x){
  mean(x[is.finite(x)], na.rm = TRUE)
}

finite.median <- function(x){
  median(x[is.finite(x)], na.rm = TRUE)
}

finite.mad <- function(x){
  mad(x[is.finite(x)], center = finite.median(x), na.rm = TRUE)
}

# Compute statistics describing data quality
# RSD: relative standard deviation: sd/mean of QCs
# RSD_r: non-parametric, robust version of RSD: MAD/median of QCs
# D_ratio: sd of QCs / sd of biological samples
# D_ratio_r: MAD of QCs / MAD of biological samples
# Detection rate: proportion of non-missing values in QC samples
# group.col should have value "QC" for QC samples
compute.quality.measures <- function(data, responses, group.col){
  require(foreach)

  qc.data <- data[data[group.col] == "QC", responses]
  sample.data <- data[data[group.col] != "QC", responses]
  results <- foreach(i=1:length(responses), .combine=rbind, .export = c("finite.sd", "finite.mad", "finite.mean", "finite.median")) %dopar% {
    response <- responses[i]
    qc.response <- qc.data[,response]
    RSD <- finite.sd(qc.response)/abs(finite.mean(qc.response))
    RSD_r <- finite.mad(qc.response)/abs(finite.median(qc.response))
    D_ratio <- finite.sd(qc.response)/finite.sd(sample.data[,response])
    D_ratio_r <- finite.mad(qc.response)/finite.mad(sample.data[,response])
    Detection_rate <- sum(!is.na(qc.response != 1))/length(qc.response)
    data.frame("Compound" = response, RSD, RSD_r, D_ratio, D_ratio_r, Detection_rate, stringsAsFactors = FALSE)
  }
  results
}

# Compute summary statistics of quality measures at each stage
# QM is a data frame with columns from quality measurements, and a column "Stage", a factor or a character indicating the stage
quality.measure.summary <- function(qm, rsd.cutoff = 0.2, d.ratio.cutoff = 0.5, detection.rate.cutoff = 0.7){
  # Compute mean and median
  smry <- qm %>%
    group_by(Stage) %>%
    summarise_at(.vars = vars(-Compound), .funs = funs(mean=finite.mean, median = finite.median))
  # Compute proportion of passes by different criteria
  smry <- cbind(smry, qm %>%
                     group_by(Stage) %>%
                     summarise(RSD_pass = sum(RSD < rsd.cutoff)/n(),
                               D_ratio_pass = sum(D_ratio < d.ratio.cutoff)/n(),
                               Detection_rate_pass = sum(Detection_rate > detection.rate.cutoff)/n(),
                               Complete_pass = sum(RSD < rsd.cutoff & D_ratio < d.ratio.cutoff & Detection_rate > detection.rate.cutoff)/n()))
  smry
}


save.quality.measure.summary.plot <- function(qm.full, output.file){
  # Set the factor levels correctly
  qm.full$Stage <- factor(qm.full$Stage, levels <- c("Original", "Drift_corrected", "Cleaned", "Imputed"))
  qm.long <- qm.full %>%
    gather("Measure", "Value", -Stage, - Compound)

  ylim1 = boxplot.stats(qm.long$Value)$stats[c(1, 5)]

  p <- ggplot(qm.long, aes(x = Measure, y = Value)) +
    geom_boxplot(aes(color = Stage)) +
    coord_cartesian(ylim = ylim1*1.05) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12)) +
    scale_color_brewer(type = "seq", palette = "Dark2") +
    labs(x = "", y = "")

  pdf(output.file, width = 12, height = 12)
  plot(p)
  dev.off()
}


# Perform drift correction from batchCorr package by Carl Brunius
# group.col: the name of the column, which should tell if a sample is a QC sample (group.col == "QC") or not
# inj.col: the name of the column which should include the injection order
# G: the possible numbers of clusters to consider
# model.names: the names of models to consider in the clustering, NULL makes the algorithm use defaults of the mclust package
# clean.responses: boolean indicating if the features should be cleaned after drift correction (recommended)
# cv.limit: the cutoff for CV in cleaning up features
# report: boolean indicating if the output plots of feature clustering should be saved
# report.dir: The directory where the clustering plots will be saved NOTE: the function will set this directory as your working directory!
correct.drift <- function(data, responses, group.col, injection.col, G = seq(1,52,by=3), model.names = NULL, clean.responses = TRUE, cv.limit = 0.2, log.file, logging, report = FALSE, report.dir = NULL){
  require(batchCorr)
  require(dplyr)

  if(!sum(data[,group.col] == "QC")){
    stop("The group column should tell if samples are QC samples or not.")
  }
  if(length(unique(data[,injection.col])) != nrow(data)){
    stop("The values in the injection number column should be unique")
  }

  # Scaling responses by sd
  scaled.data <- data
  scaled.data[responses] <- scale(scaled.data[,responses], center = FALSE, scale = apply(data[,responses], 2, sd, na.rm = TRUE))
  # Generate QC object
  qc.feats <- scaled.data[scaled.data[,group.col] == "QC", responses] # raw feature matrix
  qc.inj <- data[data[,group.col] == "QC",injection.col] # injection numbers
  qc.object <- list(inj = qc.inj,
                    Feats = qc.feats)
  # Generate object of data to be corrected (complete batch data)
  corr.feats <- scaled.data[responses]# scaled feature matrix
  corr.inj <- data[,injection.col] # injection numbers
  corr.object <- list(inj = qc.inj,
                      Feats = corr.feats)

  # Perform drift correction
  # Specify folder for pdf reports
  if(report){
    setwd(report.dir)
  }
  drift.corr.object <- clust(QCInjs = qc.object$inj, QCFeats = qc.object$Feats, report = report, modelNames = model.names, G = G) %>% # cluster responses based on drift pattern
    driftCalc(report = report, spar = NULL) %>% # compute drift functions
    driftCorr(refType = "none", CorrObj = corr.object, report = report) # Correct drift in responses where correction moves QCs closer to each other
  # Clean responses based on CV on QC samples
  if(clean.responses){
    drift.corr.object <- drift.corr.object %>%
      cleanVar(CVlimit = cv.limit, report = report)
    drift.corrected.data <- drift.corr.object$TestFeatsFinal
  } else{
    drift.corrected.data <- drift.corr.object$TestFeatsCorr
  }
  # Report numbers of responses before and after analysis
  n.resp.before <- length(responses)
  n.resp.after <- ncol(drift.corrected.data)
  n.resp.delta <- n.resp.before - n.resp.after
  log.text(paste("Responses before drift correction:", n.resp.before, "after drift correction:", n.resp.after, "removed:", n.resp.delta), log.file, logging)
  new.responses <- colnames(drift.corrected.data)
  drift.corrected.data <- as.data.frame(drift.corrected.data)
  # Scale back to original scale
  drift.corrected.data[new.responses] <- scale(drift.corrected.data[new.responses], center = FALSE, scale = apply(data[,responses], 2, function(x){1/sd(x, na.rm = TRUE)}))
  # Join drift corrected responses to data (replace original responses data with drift corrected data)
  sample.info <- data[colnames(data)[!colnames(data) %in% responses]]
  drift.corrected.data <- cbind(sample.info, drift.corrected.data)
  return(list(data = drift.corrected.data,
              responses = new.responses,
              n.responses = data.frame(before = n.resp.before,
                              after = n.resp.after,
                              delta = n.resp.delta)))
}

# Correct within batch drift by cubic spline correction
# Fits a cubic spline to QC data, then corrects the whole data for that trend
# group.col = column name, should contain QC information
# order.col = column name, should contain run order
# smooth.par: set the smoothing parameter. If NULL, will be chosen by cross validation
# smooth.par.low & high: limits for the smoothing parameter in cross validation

# correct.criterion & keep.criterion: characters giving criteria for correcting and keeping a response
# Should be compatible with dplyr::filter
# correct.criterion is applied to the change in quality measures
# keep.criterion is applied to the final quality measures.
# If the correction criterion is met, this means the post-correction quality measures. Otherwise, the original ones.
# See examples from the default ones

# Optional plotting:
# plot path: partial filename, kept and removed responses' drift correction plots will be plotted: plot.path_kept.pdf and plot.path_removed.pdf
correct.drift.spline <- function(data, responses, group.col, order.col, smooth.par = NULL, smooth.par.low = 0.4, smooth.par.high = 1.5,
                                 correct.criterion = "RSD_r < 0 & D_ratio_r < 0",
                                 keep.criterion = "(RSD_r < 0.2 & D_ratio_r < 0.4) | (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)",
                                 plotting = TRUE, plot.path = NULL, log.file = NULL, logging = NULL){
  # Set plotting file and load plotting packages
  if(plotting){
    if(any(is.null(plot.path), is.null(log.file), is.null(logging)))
      stop("Please specify plotting file, log file and logging!")
    require(ggplot2)
    require(gridExtra)
  }

  # Get QC data
  qc.data <- data[data[,group.col] == "QC",]
  qc.order <- qc.data[,order.col]
  # Overall injection order
  full.order <- data[, order.col]

  # Intializing vectors telling if the responses are corrected and/or kept
  corrected <- rep(TRUE, length(responses))
  kept <- rep(TRUE, length(responses))
  raised <- rep(FALSE, length(responses))
  # Initializing list for plots
  plots <- list()
  for(i in 1:length(responses)){
    response <- responses[i]
    # Get indexes where the response is detected in the QCs
    qc.detected <- which(!is.na(qc.data[,response]))
    if(length(qc.detected) < 0.5*length(qc.order)){
      corrected[i] <- FALSE
      next
    }
    # Fit cubic spline based on QC data
    spline.fit <- smooth.spline(x = qc.order[qc.detected], y = qc.data[qc.detected, response], spar = smooth.par, all.knots = TRUE, control.spar = list("low" = smooth.par.low, "high" = smooth.par.high))
    # Predict values for all samples
    predicted.values <- predict(spline.fit, full.order)$y
    # Sometimes predicted values can slip to the negative side. This is obviously bad, and needs to be avoided
    # The initial values are raised by their mean until the predicted values are all on the positive side.
    # NOTE that this changes the output values. Right now, we are doing a rank normalization afterwards, so it doesn't matter. But in the future, it might.
    raise <- 0
    while(any(predicted.values < 0)){
      raised[i] <- TRUE
      # Double the scaling factor
      raise <- raise + mean(qc.data[qc.detected, response])
      # Fit cubic spline based on QC data
      spline.y <- qc.data[qc.detected, response] + raise
      spline.fit <- smooth.spline(x = qc.order[qc.detected], y = spline.y , spar = smooth.par, all.knots = TRUE, control.spar = list("low" = smooth.par.low, "high" = smooth.par.high))
      # Predict values for all samples
      predicted.values <- predict(spline.fit, full.order)$y
    }
    # Compute correction factors
    correction.factors <- predicted.values[1]/predicted.values
    # Correct the drift
    corrected.response.values <- (data[,response] + raise) * correction.factors
    # Combine values into data frame
    df <- data.frame(Original = data[,response], Predicted = predicted.values - raise, Corrected = corrected.response.values)
    df[group.col] <- data[,group.col]
    df[order.col] <- data[,order.col]
    # Compute quality before and after
    qm.tmp <- compute.quality.measures(df, c("Original", "Corrected"), group.col)
    # Compute changes in quality metrics, and apply the correction criterion
    qm.diff <- data.frame(qm.tmp[2,2:5] - qm.tmp[1,2:5])
    # Apply the correction criterion
    corrected[i] <- paste0("qm.diff %>% dplyr::filter(", correct.criterion, ") %>% nrow() %>% as.logical()") %>% parse(text = .) %>% eval()

    if(plotting){
      # Plots with quality measurements
      # Draw plot of original values and fitted line, including smoothing parameter and quality measurements
      p1 <- ggplot(df, aes_string(x = order.col)) +
        geom_point(aes_string(y = "Original", color = group.col, shape = group.col)) +
        geom_line(aes(y = Predicted), color = "grey") +
        labs(title = response,
             subtitle = paste(paste(colnames(qm.tmp[2:5]), round(qm.tmp[1,2:5],digits = 3) , sep = ": ", collapse = ", "), ", smoothing parameter: ", round(spline.fit$spar, digits = 3))) +
        theme_bw()
      # Draw plot of corrected values, record if the drift correction was done or not
      p2 <- ggplot(df, aes_string(x = order.col)) +
        geom_point(aes_string(y = "Corrected", color = group.col, shape = group.col)) +
        labs(subtitle = paste0(paste(colnames(qm.tmp[2:5]), round(qm.tmp[2,2:5],digits = 3) , sep = ": ", collapse = ", "),
                              ", Drift corrected: ", corrected[i])) +
        theme_bw()
      # Combine plots vertically
      p <- arrangeGrob(p1, p2, nrow = 2)
      plots <- c(plots, list(p))
    }
    # Take either the drift corrected or original values
    # Aply keep.criterion accordingly
    if(corrected[i]){
      data[response] <- df$Corrected
      kept[i] <- paste0("qm.tmp %>% dplyr::filter( Compound == 'Corrected', ", keep.criterion, ") %>% nrow() %>% as.logical()") %>% parse(text = .) %>% eval()
    }else{
      kept[i] <- paste0("qm.tmp %>% dplyr::filter( Compound == 'Original', ", keep.criterion, ") %>% nrow() %>% as.logical()") %>% parse(text = .) %>% eval()
    }
  }
  # Plot the kept and removed responses' drift correction plots to separate files
  if(plotting){
    pdf(plot.path, onefile = TRUE)
    for(j in 1:length(plots)){
      plot(plots[[j]])
    }
    dev.off()
    log.text(paste("Saved drift correction plots to:", plot.path), log.file = log.file, logging = logging)
  }
  # Return data with drift corrected values and a data frame describing which responses were corrected and/or kept
  return(list(drift.corrected.data = data, operations = data.frame(response = responses, corrected, kept, raised, stringsAsFactors = FALSE)))
}

save.drift.correction.plots <- function(data, responses, group.col, order.col, smooth.par = NULL, smooth.par.low = 0.5, smooth.par.high = 1.5, plot.file, width = 10, height = 10, log.file, logging){


  # Get QC data and QC injection order
  qc.data <- data[data[,group.col] == "QC",]
  qc.order <- qc.data[,order.col]
  # Overall injection order
  full.order <- data[, order.col]

  plots <- foreach(i = 1:length(responses), .combine = c, .packages = c("ggplot2", "gridExtra"), .export = c("finite.sd", "finite.mad", "finite.mean", "finite.median", "compute.quality.measures")) %do% {
    response <- responses[i]
    # Get indexes where the response is detected in the QCs
    qc.detected <- which(!is.na(qc.data[,response]))
    # Fit cubic spline based on QC data
    spline.fit <- smooth.spline(x = qc.order[qc.detected], y = qc.data[qc.detected, response], spar = smooth.par, all.knots = TRUE, control.spar = list("low" = smooth.par.low, "high" = smooth.par.high))
    # Predict values for all samples
    predicted.values <- predict(spline.fit, full.order)$y
    # Compute correction factors
    correction.factors <- predicted.values[1]/predicted.values
    # Correct the drift
    corrected.response.values <- data[,response] * correction.factors
    # Gather results
    df <- data.frame(Original = data[,response], Predicted = predicted.values, Corrected = corrected.response.values)
    df[group.col] <- data[,group.col]
    df[order.col] <- data[,order.col]
    # Compute quality before and after
    qm.tmp <- compute.quality.measures(df, c("Original", "Corrected"), group.col)
    # Draw plot of original values and fitted line
    p1 <- ggplot(df, aes_string(x = order.col)) +
      geom_point(aes_string(y = "Original", color = group.col, shape = group.col)) +
      geom_line(aes(y = Predicted), color = "grey") +
      labs(title = response,
           subtitle = paste(paste(colnames(qm.tmp[2:5]), round(qm.tmp[1,2:5],digits = 3) , sep = ": ", collapse = ", "), ", smoothing parameter: ", round(spline.fit$spar, digits = 3))) +
      theme_bw()
    # Draw plot of corrected values
    p2 <- ggplot(df, aes_string(x = order.col)) +
      geom_point(aes_string(y = "Corrected", color = group.col, shape = group.col)) +
      labs(subtitle = paste(colnames(qm.tmp[2:5]), round(qm.tmp[2,2:5],digits = 3) , sep = ": ", collapse = ", ")) +
      theme_bw()
    # Combine plots vertically
    p <- arrangeGrob(p1, p2, nrow = 2)
    list(p)
  }

  pdf(plot.file, width = width, height = height)
  for(p in plots){
    plot(p)
  }
  dev.off()
  log.text(paste("Saved drift correction plots to:", plot.file), log.file = log.file, logging = logging)
}

# Impute missing values from data with random forest
# values of 1 interpreted as missing values
# wrapping in memoise is useful for faster testing
impute.random.forest <- memoise(function(data, responses, parallelize = FALSE, log.file, logging){
  require(missForest)

  data.responses <- data[responses]
  n.ones <- sum(data.responses == 1, na.rm = TRUE)
  if(n.ones){
    log.text(paste("Found", n.ones, "values of 1 in the data, setting them to NA"))
    data.responses[data.responses == 1] <- NA
  }
  parallelize <- ifelse(parallelize, "variables", "no")

  mf <- missForest(xmis = data.responses, parallelize = parallelize)
  imputed.responses <- mf$ximp
  log.text(paste("\n Out-of-bag error in random forest imputation:", round(mf$OOBerror, digits = 3), "\n") , log.file=log.file, logging=logging)
  data[responses] <- imputed.responses
  data
})

# Find out if the signals are present in at least 80% (adjustable) of samples in any non-QC group
detected.in.groups <- function(data, responses, id.col, group.col, detection.limit = 0.8){
  df <- data[, c(id.col, group.col,responses)] %>%
    tidyr::gather_("Response", "Intensity", responses) %>%
    dplyr::group_by_("Response", group.col) %>%
    dplyr::summarise(proportion_found = sum(!is.na(Intensity))/n()) %>%
    tidyr::spread(group.col, "proportion_found")
  df.proportions <- df[-1]
  df.proportions["QC"] <- NULL
  df$Detected <- apply(df.proportions, 1, function(x){any(x > detection.limit)})
  df
}

to.original.format <- function(x){
  require(dplyr)
  x %>%
    gsub("COMPOUND_", "", .) %>%
    gsub("_", ".", .) %>%
    gsub("a", "@", .)
}

# Convert COMPOUND_123_456a3_456 to original format 123.456@3.456
convert.compound.format <- function(x){
  sapply(x, to.original.format) %>% unname()
}


check_imputation <- function(orig.data, imputed.data, responses){
  # Flatten the data matrices into vectors for easy comparison
  orig <- orig.data[responses] %>% as.matrix() %>% as.vector()
  imputed <- imputed.data[responses] %>% as.matrix() %>% as.vector()
  
  check <- is.na(orig) | orig == imputed
  data.frame(response = responses, check = check)
}

