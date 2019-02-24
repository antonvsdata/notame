library(nlme)
library(ggplot2)

##################
# Visualizations #
##################

##########################################################
# 0. New functions with with 1 and 2 combined            #
# 1. save.x.plot: prepare data, compute, draw/plot, save #
# 2. compute, draw/plot functions                        #
##########################################################


# Parse new vector (column) from available choices.
## df: data.frame
## id.col.name: string
## id.col.abbrev: string or null
## group.names: list or null
## group.abbrevs: list or null
parse.new.id.col <- function(df, id.col.name, id.col.abbrev=NULL, group.names=NULL, group.abbrevs=NULL){
  
  new.id.col <- NULL
  # group names are not available
  if(is.null(group.names)){
    # id column abbreviation is not available
    if(is.null(id.col.abbrev)){
      new.id.col <- paste(id.col.name, df[,id.col.name], sep="")
      # id column abbreviation is available
    } else {
      new.id.col <- paste(id.col.abbrev, df[,id.col.name], sep="")
    }
    # group abbreviations not available, group names are available
  } else if(is.null(group.abbrevs) & !is.null(group.names)){
    # id column abbreviation is not available
    if(is.null(id.col.abbrev)){
      for(g in group.names){
        new.id.col <- paste(new.id.col, g, df[,g], sep="")
      }
      # add id column original dataframe name
      new.id.col <- paste(id.col.name, df[,id.col.name], new.id.col, sep="")
      # use id column abbreviation
    } else {
      for(g in group.names){
        new.id.col <- paste(new.id.col, g, df[,g], sep="")
      }
      # add id column original dataframe name
      new.id.col <- paste(id.col.abbrev, df[,id.col.name], new.id.col, sep="")
    }
  # group abbreviations are available
  } else if(!is.null(group.abbrevs) & !is.null(group.names)){
    # check group names and abbrevs are the same length
    if(length(group.names) == length(group.abbrevs)){
      # id column abbreviation is not available
      if(is.null(id.col.abbrev)){
        for(g in 1:length(group.names)){
          new.id.col <- paste(new.id.col, group.abbrevs[g], df[,group.names[g]], sep="")
        }
        new.id.col <- paste(id.col.name, df[,id.col.name], new.id.col, sep="")
      # id column abbreviation is available
      } else{
        for(g in 1:length(group.names)){
          new.id.col <- paste(new.id.col, group.abbrevs[g], df[,group.names[g]], sep="")
        }
        new.id.col <- paste(id.col.abbrev, df[,id.col.name], new.id.col, sep="")
      }
      
    } else{
      stop("Mismatch between the numbers of group.names and group.abbrevs parametres! Make sure group.names and -abbrevs are lists of same length!")
    }
  }
  # return new id column (vector)
  new.id.col
}


# Create start, end intervals according to page.limit and lenght of given vector and output.filename_x
## eg. 7501 pages with 1000 intervals(page.limit) would start from 1 to (end)1000, next: (start)1001 to (end)2000... (start)7000 to (end)7501
divided.plotting <- function(page.limit, responses, output.filename){
  
  listlength <- length(responses)
    
  # Paging per document .pdf
  tmp.df <- data.frame(n=seq(1, ceiling(listlength/page.limit), 1))
  tmp.df$start <- c(1, seq(page.limit+1, listlength, page.limit))
  if(listlength == tail(seq(page.limit, listlength, page.limit), n=1)){
    tmp.df$end <- seq(page.limit, listlength, page.limit)
  } else {
    tmp.df$end <- c(seq(page.limit, listlength, page.limit), listlength)  
  }
  # output file names: output_XofY.pdf
  tmp.df$output <- paste0(output.filename,tmp.df$n, "of", length(tmp.df$n), ".pdf")

  tmp.df
}

# tSNE
save.tsne.plot <- function(data, responses, group.col, output.file, logging, log.file){
  require(Rtsne)
  
  tmp.matrix <- as.matrix(data[responses])
  set.seed(42)
  data.rtsne <- Rtsne(tmp.matrix)
  
  #plot(dat.rtsne$Y, col=data$GROUP)
  
  # Dataframe for ggplot
  temp.df <- data.frame(data.rtsne$Y, GROUP=data[,group.col])
  # Change name of the (grouping) column
  colnames(temp.df) <- c("X1","X2",group.col)
  
  # Plot with ggplot2
  p <- ggplot(temp.df, aes_string(x="X1", y="X2", colour=group.col, shape=group.col, group=group.col)) +
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

######
# 0. #
######

corr.heatmap <- function(data, responses, phenotypes, corr.file.prefix, responses.title=NULL, phenotypes.title=NULL, rm.axis.labels=FALSE, plot.title="Heatmap - Correlation", output.file=NULL, logging, log.file, grad.colours=NULL, grad.breaks=NULL, grad.limits=NULL, pdf.height.divisor=1, pdf.width.divisor=1, distance.method="euclidean", hclust.method="complete") {
  require(reshape)
  require(gridExtra)
  require(ggplot2)
  
  tempdata <- data.frame()
  
  for (phenotype in phenotypes ) {
    # generate data file name
    pheno.file <- paste(corr.file.prefix, phenotype, "_correlation.csv", sep="")
    # read data from file (saved correlation results)
    tryCatch({
      read_corr_data <- read.csv(pheno.file, header=TRUE)
    }, error = function(e) {
      print(e$message)
      log.text(paste("Error, could not open the following file: ", pheno.file, "\nMake sure correlation files for data are created.", sep=""), log.file=log.file, logging=logging)
    })
    # Add phenotype information to data.frame
    read_corr_data$Pheno <- phenotype
    # Subset only relevant data
    read_corr_data <- data.frame(read_corr_data$Pheno, read_corr_data$Trait, read_corr_data$Estimate, read_corr_data$p.value_adjusted_BH)

#     # cutoff by p-value if given
#     if(!is.null(p.cutoff)){
#       read_corr_data <- read_corr_data[read_corr_data[,"read_corr_data.p.value_adjusted_BH"] < p.cutoff, c("read_corr_data.Pheno", "read_corr_data.Trait", "read_corr_data.Estimate")]
#     }

    # transform/transpose data, phenotype is on single row and responses as columns
    pheno.row <- cast(read_corr_data, read_corr_data.Pheno~read_corr_data.Trait, value="read_corr_data.Estimate")
    
    # bind row to existing data.frame
    tempdata <- rbind(tempdata, pheno.row)
  }  

  corr.data <- cor(na.omit(data))

#   print("f0lumbmd, Gp")
#   print(corr.data["f0lumbmd", "Gp"])
#   print("cor.test")
#   print(cor.test(data[,"f0lumbmd"], data[,"Gp"]))
#   print("tempdata")
#   rownames(tempdata) <- tempdata$read_corr_data.Pheno
#   print(tempdata["f0lumbmd","Gp"])
#  print(head(tempdata))

  filu <- paste(corr.file.prefix, "_correlation_TEST.csv", sep="")
  write.csv(corr.data, file=filu, row.names=TRUE, na="")
  

  # Subset data to group 1 (phenotypes) and group 2 (responses) for each independent distance and hclust, rows and columns
  pheno.corr.data <- corr.data[phenotypes, phenotypes]
  resp.corr.data <- corr.data[responses, responses]
  
  # Hiearchical clustering, distances
  pheno.dist <- dist(pheno.corr.data, method=distance.method)
  resp.dist <- dist(resp.corr.data, method=distance.method)
  
  # hiearchical clustering
  hc.pheno <- hclust(pheno.dist, method=hclust.method)
  hc.resp <- hclust(resp.dist, method=hclust.method)

  # extracting row and col order from hiearchical clustering
  hc.pheno.order <- hc.pheno$labels[c(hc.pheno$order)]
  hc.resp.order <- hc.resp$labels[c(hc.resp$order)]

  # Unique rows (for phenotypes) and cols (for responses)
  unequal.data <- corr.data[phenotypes, responses]

#   print(hc.pheno.order)
#   print(hc.resp.order)
  
  #melted <- melt(tempdata)
  
  melted <- melt(unequal.data)

  ## if both column names / titles are given NULL, then remove axis labels from plot
  if(is.null(phenotypes.title) && is.null(responses.title)){ rm.axis.labels <- TRUE }
  # data column names are also used as titles, they need to be defined even if given NULL
  if(is.null(responses.title)){ responses.title <- "responses" }
  if(is.null(phenotypes.title)){ phenotypes.title <- "phenotypes" }
  
  colnames(melted) <- c(phenotypes.title, responses.title, "Correlation")

  melted[,phenotypes.title] <- factor(melted[,phenotypes.title], levels=hc.pheno.order, ordered=TRUE)
  melted[,responses.title] <- factor(melted[,responses.title], levels=hc.resp.order, ordered=TRUE)

  p <- draw.heatmap(melted, responses.title, phenotypes.title, fill.column="Correlation", title=plot.title, remove.axis.labels=rm.axis.labels, x.axis.text.angle=90, value.limits=c(1,-1), md=0, gradient.colours=grad.colours, gradient.breaks=grad.breaks, gradient.limits=grad.limits)

  # Save plot to file
  if(!is.null(output.file)){
    # Approximate height and width
    approx.height <- length(unique(melted[,phenotypes.title]))
    approx.width <- length(unique(melted[,responses.title]))
    
#     print(approx.height)
#     print(approx.width)
#     print(round(approx.width/pdf.width.divisor))
#     print(round(approx.height/pdf.height.divisor))

    # Set width and height for saving pdf file
    ## difficult to automate, pdf.width and height.divisor are given as parametres
    pdf(output.file, width=round(approx.width/pdf.width.divisor), height=round(approx.height/pdf.height.divisor))
    # Print to arrangeGrob to device
    grid.arrange(p)
    # Turn of device
    dev.off()
    
    if(!is.null(logging) & !is.null(log.file)){
      # Log 
      log.text(paste("Saved correlation heatmap phenotypes vs responses - plot to: ", output.file, sep=""), log.file=log.file, logging=logging)
    } else{
      stop("Error, logging and log.file need to be given in order to save to file!")
    }
  } else{
    # Return plot (no output.file defined)
    p
    if(!is.null(logging) & !is.null(log.file)){
      # Log 
      log.text(paste("Produced heatmap phenotypes vs responses - plot to variable"), log.file=log.file, logging=logging)
    }
  }
  # return plot even if no output.file defined
  p
}

# Changes?: x.col to list of variables (but max 2.?; how to visualize 
line.plots <- function(data, responses, x.col, x.levels=NULL, group.col=NULL, group.col.title=NULL, group.levels=NULL, limit.amount=NULL, position.dodge=0.3, log.file=NULL, logging=NULL, output.file=NULL) {
  
  require(doMC)
  require(ggplot2)
  
  tmp.data <- data

  if(!is.null(x.levels)){
    tmp.data <- droplevels(tmp.data[tmp.data[,x.col] %in% x.levels,])
  } else { # is NULL
    x.levels <- levels(data[,x.col])
  }
  
  if(!is.null(group.col)){
    # Select only data that matches group.levels in column group.col
    if(!is.null(group.levels)){
      tmp.data <- droplevels(tmp.data[tmp.data[,group.col] %in% group.levels,])
    } else { # is NULL
      group.levels <- levels(data[,group.col])
    }
  }
  
  # Plotting
  
  # Mean point with errorbars, multiple groups with dodge 
  pd <- position_dodge(position.dodge)
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  plots <- list()
  
  if(is.null(limit.amount)){
    limit.amount <- length(responses)
  }
  
  # group.col to group.col.title if group.col.title is not given
  if(!is.null(group.col)){
    if(is.null(group.col.title)){
      group.col.title <- group.col
    }
  }
  
  plots <- foreach(i=1:limit.amount, .combine=c) %dopar% {
    library(ggplot2)
    
    response <- responses[i]
    print(response)
    print(i)
    
    if(!is.null(group.col)){  
      p <- ggplot(tmp.data, aes_string(x=x.col, y=response, group=group.col, colour=group.col)) +
        # Errorbars with solid lines
        stat_summary(fun.data = "mean_se",
                     geom="errorbar",position=pd, width=.5,
                     fun.ymin = function(x) mean(x) - sd(x), 
                     fun.ymax = function(x) mean(x) + sd(x)) +
        # Plot point to mean
        stat_summary(fun.y = mean,
                     geom = "point",
                     position=pd, size=3,
                     fun.ymin = min,
                     fun.ymax = function(x) mean(x) + sd(x),
                     aes_string(shape=group.col)) +
        # Line from mean to mean between for example timepoints
        stat_summary(fun.y = mean,
                     geom = "line",
                     position=pd, size=0.5,
                     fun.ymin = min,
                     fun.ymax = function(x) mean(x) + sd(x),
                     aes_string(linetype=group.col))
      
      p <- p + scale_colour_discrete(name=group.col.title) +
        scale_shape_discrete(name=group.col.title) +
        scale_linetype_discrete(name=group.col.title) +
        scale_fill_discrete(name=group.col.title)
    
    # No group.col input
    } else {
      p <- ggplot(tmp.data, aes_string(x=x.col, y=response)) +
        # Errorbars with solid lines
        stat_summary(fun.data = "mean_se",
                     geom="errorbar",position=pd, width=.5,
                     fun.ymin = function(x) mean(x) - sd(x), 
                     fun.ymax = function(x) mean(x) + sd(x)) +
        # Plot point to mean
        stat_summary(fun.y = mean,
                     geom = "point",
                     position=pd, size=3,
                     fun.ymin = min,
                     fun.ymax = function(x) mean(x) + sd(x)
                     ) +
        # Line from mean to mean between for example timepoints
        stat_summary(fun.y = mean,
                     geom = "line",
                     position=pd, size=0.5,
                     fun.ymin = min,
                     fun.ymax = function(x) mean(x) + sd(x),
                     aes_string(group=1),
                     linetype=2)
    }
    
    p <- p + theme(legend.position="top")
    p <- p + ggtitle(response) + 
      theme(plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle=90, hjust=0.5))
    p <- p + theme_bw()
    
    p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)
    
    list(p)
  }
  
  if(!is.null(output.file)){
    if(is.null(log.file) | is.null(logging)){ stop("Error, log.file and logging required for saving to file!") }
    # Save plotlist and log
    save.plotlist.pdf(output.file, plots, log.description="Saved line plots to: ", log.file, logging)
    
  } else{
    plots
    
    if(!is.null(logging) & !is.null(log.file)){
      # Log 
      log.text(paste("Produced list of lineplots - plot to variable"), log.file=log.file, logging=logging)
    }
  }
  plots
}

# Response heatmap
## Extract order of metabolites from metabolite subset by hierarchical clustering,
## create a combined group variable from Gest_group and Timepoint,
## aggregate data by function (default is mean,)
## draw heatmap from aggregated data
##
## colours: high: [1], mid: [2], low: [3]
heatmap.response <- function(data, responses, group.names, group.abbrevs=NULL, x.label=NULL, y.label="Response", fun="mean", colours=c("red","orange","white"), x.text.angle=0, log.file=NULL, logging=NULL, output.file=NULL, width=NULL, height=NULL, scale=0.3){

  require(reshape)
  require(gridExtra)
  require(ggplot2)
  
  tmp.data <- data
  
  # Create new.id.column based on group.names (and group.abbrevs if defined)
  new.id.col <- NULL
  
  # No group.abbrevs defined
  if(is.null(group.abbrevs)){
    for(group in group.names){
      new.id.col <- paste(new.id.col, group, tmp.data[,group], sep="")
    }
  # use group abbreviations
  } else {
    if(length(group.names) == length(group.abbrevs)){
      for(g in 1:length(group.names)){
        group <- group.names[g]
        new.id.col <- paste(new.id.col, group.abbrevs[g], tmp.data[,group], sep="")
      }
    } else{
      stop("Error, input parametres group.names and group.abbrevs not same length!")
    }
  }
  
  # make it a factor
  tmp.data$combined <- as.factor(new.id.col)
  tmp.data <- na.omit(tmp.data[c("combined", group.names, responses)])
  tmp.data <- droplevels(tmp.data)
  
  data.dist <- dist(t(tmp.data[responses]))
  #data.dist <- dist(data[responses])
  hc <- hclust(data.dist)
  hc.order <- hc$labels[c(hc$order)]
  
  # Generate x.label if not given
  if(is.null(x.label)){
    for(group in group.names){
      x.label <- paste(x.label, group, sep="")
      x.label <- paste(x.label, "_", sep="")
    }
  }
  
  # Aggregate data by mean of the metabolite per combined group and time factor (GT), with function mean as default
  melted <- melt.responses.by.groups(tmp.data, responses, "combined", x.label=x.label, y.label=y.label, fun=fun)
  # reorder Metabolites to hclust order
  melted$Response <- factor(melted$Response, levels=hc.order, ordered=TRUE)  
  # Plot
  hm.plot <- draw.heatmap(melted, x.label=x.label, y.label=y.label, x.axis.text.angle=x.text.angle, colour.low=colours[3], colour.mid=colours[2], colour.high=colours[1])
  
  # Save to output.file
  if(!is.null(output.file)){
    if(is.null(logging) & is.null(log.file)){
      stop("Error, logging and log.file need to be given inorder to save to file!")
    } else{
      if(is.null(width) & is.null(height)){
        h <- convertHeight(grobHeight(tableGrob(hm.plot)), "in", valueOnly=TRUE)
        print(h)
        w <- convertWidth(grobWidth(tableGrob(hm.plot)), "in", valueOnly=TRUE)
        print(w)
        ggsave(hm.plot, file=output.file, scale=scale, width=w, height=h, limitsize=FALSE)
      } else if(width == FALSE & height == FALSE){
        ggsave(hm.plot, file=output.file, limitsize=FALSE)
      } else {
        ggsave(hm.plot, file=output.file, limitsize=FALSE, height=height, width=width)
      }
      graphics.off()
      log.text(paste("Saved subset heatmap metabolite grouping by group and time heatmap - plot to: ", output.file, sep=""), log.file=log.file, logging=logging)
    }
  }
  hm.plot
}

# Hiearchical clustering samples vs samples (distance: default: euclidean)
heatmap.samples <- function(data, responses, id.col.name, x.label="", y.label="", group.bar.label="", log.file=NULL, logging=NULL, id.col.abbrev=NULL, separate.group.bar=NULL, group.names=NULL, group.abbrevs=NULL, output.file=NULL, dist.method="euclidean", colours=c("white", "orange", "red")){

  require(reshape)
  require(gridExtra)
  require(ggplot2)
  
  # Prepare data.frame
  tmp.data <- data  
  ## Create new id column based on available input parametres
  tmp.data$new.id <- parse.new.id.col(df=tmp.data, id.col.name=id.col.name, id.col.abbrev=id.col.abbrev, group.names=group.names, group.abbrevs=group.abbrevs)
  ## omit naslog
  tmp.data <- na.omit(tmp.data[c("new.id", id.col.name, group.names, responses)]) # TEST!

  # Set row names to new.id
  rownames(tmp.data) <- tmp.data$new.id
  # Compute sample distances
  sample.distances <- dist(tmp.data[responses], method=dist.method)

  # Get order by hiearchical clustering 
  hc <- hclust(sample.distances)
  hc.order <- hc$labels[c(hc$order)]
  
  distances <- as.matrix(sample.distances)
  melted <- melt(distances)
  
  # Merge with row information and by defined group (separate.group.bar) for additional group bar
  if(is.null(separate.group.bar)){
    melted<- merge(melted, tmp.data[c("new.id")], by.x="X1", by.y="new.id")  
  } else {
    melted<- merge(melted, tmp.data[c("new.id", separate.group.bar)], by.x="X1", by.y="new.id")  
  }
  
  # reorder Samples to hclust order
  melted$X1 <- factor(melted$X1, levels=hc.order, ordered=TRUE)
  melted$X2 <- factor(melted$X2, levels=hc.order, ordered=TRUE)
  
  print("Starting generating heatmap plot...")
  
  p <- draw.heatmap(melted, x.axis.text.angle=90,
                    title=expression(atop("Heatmap", atop(italic("Samples ordered by hierarchical clustering. Euclidean distance between samples."), ""))),
                    colour.low=colours[3],
                    colour.mid=colours[2],
                    colour.high=colours[1])
  p <- p + scale_x_discrete(name=x.label) +
            scale_y_discrete(name=y.label)
  
  # To draw additional group bar
  if(!is.null(separate.group.bar)){  
    # draw coloured bar
    p0 <- ggplot(data=melted, aes_string(x="X1", y=1)) +
      geom_tile(aes_string(fill=separate.group.bar)) + 
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks=element_blank(),
            axis.text=element_blank())
    
    p0 <- p0 + scale_x_discrete(name=group.bar.label)
    # scale_y_discrete(limits=0.01)
    # scale_y_continuous(expand=c(0,0))
    
    gp1 <- ggplotGrob(p)
    gp2 <- ggplotGrob(p0)
    
    maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
    gp1$widths[2:5] <- as.list(maxWidth)
    gp2$widths[2:5] <- as.list(maxWidth)
    
    
    p <- arrangeGrob(gp1,gp2, ncol=1, heights=c(10/11,1/11))
  }
  
  print("Generated heatmap")
  print(p)
  print("Starting to save to file..")
  
  # Save plot to file
  if(!is.null(output.file)){
    # Approximate samples count
    approx.count <- length(tmp.data$new.id)
    #print(approx.count)
    
    # Set width and height for saving pdf file
    pdf(output.file, width=approx.count/5, height=(approx.count/5)+2)
    # Print to arrangeGrob to device
    grid.arrange(p)
    # Turn of device
    dev.off()
    
    if(!is.null(logging) & !is.null(log.file)){
      # Log 
      log.text(paste("Saved subset heatmap samples vs samples - plot to: ", output.file, sep=""), log.file=log.file, logging=logging)
    } else{
      stop("Error, logging and log.file need to be given in order to save to file!")
    }
  } else{
    # Return plot (no output.file defined)
    p
    if(!is.null(logging) & !is.null(log.file)){
      # Log 
      log.text(paste("Produced heatmap samples vs samples - plot to variable"), log.file=log.file, logging=logging)
    }
  }
  p
}

factor.boxplots <- function(data, phenotypes, phenotypes.labels=NULL, responses, p.cutoff=NULL, corr.file.prefix, output.prefix, log.file, logging){
  require(doMC)
  require(ggplot2)
  require(reshape)
  #require(scales)
  
  for(i in 1:length(phenotypes)){
  
    phenotype <- phenotypes[i]
    
    tmp.data <- na.omit(data[c(phenotype, responses)])
    
    # Set phenotype.label from phenotypes.labels OR if NULL then set phenotype.label to phenotype (variable name)
    phenotype.label <- phenotypes.labels[i]
    if(is.null(phenotype.label)){
      phenotype.label <- phenotype
    }
    # Plotting
    cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
    
#     # perform correlation
#     results <- perform.cor(data=tmp.data, responses=responses, phenotype=phenotype)
    
    # generate data file name
    pheno.file <- paste(corr.file.prefix, phenotype, "_correlation.csv", sep="")
    # read data from file (saved correlation results)
    tryCatch({
      results <- read.csv(pheno.file, header=TRUE)
    }, error = function(e) {
      print(e$message)
      log.text(paste("Error, could not open the following file: ", pheno.file, "\nMake sure correlation files for data are created.", sep=""), log.file=log.file, logging=logging)
    })

    # adjusted p-value LOWER cutoff, if cutoff given
    if(!is.null(p.cutoff)){
      p.responses <- as.vector(results$Trait[results[,"p.value_adjusted_BH"] < p.cutoff])
    }
    if(is.null(p.cutoff)){
      p.responses <- responses
    }
    
    plots <- list()
  
    if(length(p.responses) > 0){

# For DEBUGGING:    
#       print(phenotype)
#       print(p.responses)
# 
#        res <- results[results[,"p.value_adjusted_BH"] < p.cutoff,]
#        write.csv(res, file=paste(output.prefix, phenotype, "_boxplots.csv", sep=""), row.names=FALSE, na="")
      
      # Parallelized 
      plots <- foreach(i=1:length(p.responses), .combine=c) %dopar% { 
        require(ggplot2)
        require(reshape)
        
        response <- p.responses[i]
        subset.data <- tmp.data[,c(phenotype, response)]
        
        # Set to factor if only few unique data points
        if(length(unique(tmp.data[,phenotype])) < 5){
          subset.data[,phenotype] <- factor(subset.data[,phenotype])
        }
        
        # Actual plot
        p <- ggplot(subset.data, aes_string(x=phenotype, y=response, fill=phenotype))
        p <- p + geom_boxplot()
        
        # Theme, title, texts...
        p <- p + theme(legend.position="top") + xlab(phenotype.label)
        p <- p + ggtitle(response) + 
          theme(plot.title = element_text(face="bold"),
                axis.text.x = element_text(angle=90, hjust=0.5))
        p <- p + theme_bw()
        
        list(p)
      }
      
      output.file <- paste(output.prefix, phenotype, "_boxplot.pdf", sep="")
    
      # Save plotlist and log
      tryCatch({
        save.plotlist.pdf(output.file, plots, log.description="Saved boxplots to: ", log.file, logging)
      }, error = function(e) {
        print(e$message)
        log.text(paste("Error, could not save the following file: ", output.file, sep=""), log.file=log.file, logging=logging)
      })
    } else {
      log.text(paste(phenotype, " -phenotype had 0 correlations with responses under adjusted signifance cutoff ", p.cutoff, ".\n", sep=""), log.file=log.file, logging=logging)
    }
  }
}

######
# 1. #
######

# LEGACY 
save.samplesheatmap.plot <- function(data, responses, idcol.name, idcol.abbrev, group.a.name, group.a.abbrev, group.b.name, group.b.abbrev, separate.col, separate.col.label, output.file, log.file, logging){
  require(gridExtra)
  require(ggplot2)
  # Heatmap samples vs samples
  samples.heatmap.plot <- samples.heatmap(data=data,
                                          responses=responses,
                                          id.colname = idcol.name,
                                          id.col.abbrev = idcol.abbrev,
                                          group.a.colname = group.a.name,
                                          group.a.abbrev = group.a.abbrev,
                                          group.b.colname = group.b.name,
                                          group.b.abbrev = group.b.abbrev,
                                          separate.group.bar = separate.col,
                                          separate.group.bar.label = separate.col.label,
                                          colour.high="white",
                                          colour.mid="orange",
                                          colour.low="red")
  
  # Approximate samples count
  unique.ids <- paste(idcol.abbrev, data[,idcol.name],
                           group.a.abbrev, data[,group.a.name],
                           group.b.abbrev, data[,group.b.name],
                           sep="")
  approx.count <- length(unique.ids)
  #print(approx.count)
  
  # Set width and height for saving pdf file
  pdf(output.file, width=approx.count/5, height=(approx.count/5)+2)
  # Print to arrangeGrob to device
  grid.arrange(samples.heatmap.plot)
  # Turn of device
  dev.off()
  
  # Log 
  log.text(paste("Saved subset heatmap samples vs samples - plot to: ", output.file, sep=""), log.file=log.file, logging=logging)

}

save.metabolite.heatmap <- function(data, responses, group.a.name, group.a.abbrev, group.b.name, group.b.abbrev, comb.group.label, output.file, log.file, logging) {
  require(gridExtra)
  
  hm.plot <- metabolite.grouping.heatmap(data=data, responses=responses,
                                         group.a.colname=group.a.name, group.a.abbrev=group.a.abbrev,
                                         group.b.colname=group.b.name, group.b.abbrev=group.b.abbrev,
                                         combined.group.label=comb.group.label,
                                         colour.low="white", colour.mid="orange", colour.high="red",
                                         x.text.angle=45)
  hm.plot
  h <- convertHeight(grobHeight(tableGrob(hm.plot)), "in", valueOnly=TRUE)
  #print(h)
  w <- convertWidth(grobWidth(tableGrob(hm.plot)), "in", valueOnly=TRUE)
  #print(w)
  
  ggsave(hm.plot, file=output.file, scale=0.3, width=w+10, height=h, limitsize=FALSE)
  graphics.off()
  log.text(paste("Saved subset heatmap metabolite grouping by group and time heatmap - plot to: ", output.file, sep=""), log.file=log.file, logging=logging)
}


save.box.plot <- function(data, responses, idcol.name, idcol.abbrev, group.a.name, group.a.abbrev, group.b.name, group.b.abbrev, zoom.boxplot=FALSE, output.file, log.file, logging){
  # Create plot
  p <- draw.boxplots(data=data, responses=responses,
                     id.colname = idcol.name,
                     id.col.abbrev = idcol.abbrev,
                     group.a.colname = group.a.name,
                     group.a.abbrev = group.a.abbrev,
                     group.b.colname = group.b.name,
                     group.b.abbrev = group.b.abbrev,
                     zoom.boxplot=TRUE)
  # ggSave plot, logging in ggsave.plot
  ggsave.plot(p, plot.file=output.file, log.description="Saved samples boxplots - plot to: ", log.file, logging, scale=2.5)  
}

save.line.plot <- function(data, responses, group.col, group.col.title=NULL, group.levels=NULL, x.col, x.levels=NULL, position.dodge=0.3, output.file, log.file, logging){
  
  # Create plotlist
  plots <- draw.plots(data=data, responses=responses, group.col=group.col, group.col.title=group.col.title, group.levels=group.levels, x.col=x.col, x.levels=x.levels, position.dodge=position.dodge)
  # Save plotlist and log
  save.plotlist.pdf(output.file, plots, log.description="Saved line plots to: ", log.file, logging)
}

save.histogram.plot <- function(data, responses, bin.width=30, output.file, log.file, logging){
  
  # Create plotlist
  plots <- draw.histograms(data=data, responses=responses, bin.width=bin.width)
  # Save plotlist and log
  save.plotlist.pdf(output.file, plots, log.description="Saved histogram plots to: ", log.file, logging)
}

save.pca.plot <- function(data, responses, group.a.ref, group.a.title, sampleid.labels.ref, width=6, height=6, fixed_coordinates=FALSE, output.file, log.file, logging){
  require(gridExtra)
  
  pca.result <- compute.pca(data[responses])
  
  # Same PCA plot without sample IDs
  pca.plot <- plot.pca(pca.data=pca.result, first.group=group.a.ref, first.group.title=group.a.title, fixed_coordinates=fixed_coordinates)
  
  # PCA plot with designated sample IDs in plot
  if(!is.null(sampleid.labels.ref)){
    pca.plot.2 <- plot.pca(pca.data=pca.result, first.group=group.a.ref, first.group.title=group.a.title, sampleid.labels=sampleid.labels.ref)
    p <- arrangeGrob(pca.plot, pca.plot.2, nrow=2)
  } else { 
    p <- pca.plot
  }

  #ggsave.plot(p, plot.file=output.file, log.description="Saved PCA plot to: ", log.file, logging, scale=1.5, width=6, height=12)
  save.plot(p, plot.file=output.file, log.description="Saved PCA plot to: ", log.file=log.file, logging=logging, width=width, height=height)
}

save.scatter.plot <- function(data, responses, phenotypes, pheno.label, group.col=NULL, group.col.label=NULL, group.stat_smooth=FALSE, output.file, log.file, logging){
  # Create plotlist
  plots <- draw.scatterplots(data=data, responses=responses, phenotypes=phenotypes, phenotypes.labels=pheno.label, group.col=group.col, group.col.label=group.col.label, group.stat_smooth=group.stat_smooth)

  # Save plotlist and log
  tryCatch({
    save.plotlist.pdf(output.file, plots, log.description="Saved scatter plots to: ", log.file, logging)
  }, error = function(e) {
    print(e$message)
    log.text(paste("Error, could not save the following file: ", output.file, sep=""), log.file=log.file, logging=logging)
  })
}

save.pca.faceted.plot <- function(data, responses, group.a.ref, group.a.title, sampleid.labels.ref, facet.group.ref, width=23, height=12, fixed_coordinates=FALSE, output.file, log.file, logging){
  require(gridExtra)
  
  pca.result <- compute.pca(data[responses])
  # Facets by x from previous PCA for all the data
  ###################################################
  grouped.pca.plot <- plot.pca(pca.data=pca.result, first.group=group.a.ref, first.group.title=group.a.title, facet.group=facet.group.ref, fixed_coordinates=fixed_coordinates)
  grouped.pca.plot.2 <- plot.pca(pca.data=pca.result, first.group=group.a.ref, first.group.title=group.a.title, sampleid.labels=sampleid.labels.ref, facet.group=facet.group.ref, fixed_coordinates=fixed_coordinates)

  p <- arrangeGrob(grouped.pca.plot, grouped.pca.plot.2, nrow=2)
  save.plot(p, plot.file=output.file,
              log.description=paste("Saved PCA grouped by ", group.a.title, ", facets by ", colnames(facet.group.ref), " - plot to: ", sep=""), log.file, logging,
              scale=1.5, width=width, height=height)  
  
}

save.grouped.pca.plot <- function(data, responses, group.col, group.title, color.group.col, color.group.title, title="PCA", ncol=1, nrow=1, scale=FALSE, center=FALSE, fixed_coordinates=FALSE, sampleid.labels.col, output.file, log.file, logging){
  
  # PCA according to group
  pca.group.plots <- plot.multigroup.pca(
    data=data,
    responses=responses,
    group.subset=group.col, group.subset.title=group.title,
    group.coloring=color.group.col, group.coloring.title=color.group.title,
    title.prefix=title,
    scale=scale,
    center=center,
    fixed_coordinates=fixed_coordinates)
  
  # PCA with sampleids
  pca.group.plots.v2 <- plot.multigroup.pca(
    data=data,
    responses=responses,
    group.subset=group.col, group.subset.title=group.title,
    group.coloring=color.group.col, group.coloring.title=color.group.title,
    title.prefix=title,
    sampleid.labels=sampleid.labels.col,
    scale=scale,
    center=center,
    fixed_coordinates=fixed_coordinates)
  
  #pca.group.plots
  #pca.group.plots.v2

  #ncol <- length(unlist(pca.group.plots))
  #print(ncol)
  
  height <- 6*nrow
  width <- (15/nrow)*ncol
  
  #print(height)
  #print(width)
  
  p <- do.call("arrangeGrob", c(pca.group.plots, pca.group.plots.v2, nrow=nrow, ncol=ncol))
  save.plot(p, plot.file=output.file,
              log.description=paste("Saved PCA grouped by ", group.col, ", colored by ", color.group.col, " - plot to: ", sep=""), log.file, logging,
              scale=3, width=width, height=height)
}

######
# 2. #
######

draw.histograms <- function(data, responses, bin.width=30) {
  require(doMC)
  require(ggplot2)
  
  len <- length(responses)
  
  plotteja <- foreach(i=1:len, .combine=c) %dopar% {
    require(ggplot2)

    response <- responses[i]
    
    # for bin width
    df.range <- range(data[,response])
    
    # qplot because ggplot can't handle: p <- ggplot(data[,response], aes_string(x=response)) ...
    p <- qplot(data[,response], geom = 'blank') +   
      geom_line(aes(y = ..density..), stat = 'density') +  
      geom_histogram(aes(y = ..density..), alpha = 0.4, binwidth=df.range[2]/bin.width)                        

    p <- p + geom_density()
    p <- p + theme(legend.position="top",
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank()) +
      xlab("")
    p <- p + ggtitle(response) + 
      theme(plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle=90, hjust=0.5))
    p <- p + theme_bw()
    
    
    list(p)
  }
  plotteja

}


draw.boxplots <- function(data, responses, id.colname, id.col.abbrev, group.a.colname, group.a.abbrev, group.b.colname, group.b.abbrev,
                          # x.col="new.id",
                          x.levels=NULL,
                          limit.amount=NULL,
                          zoom.boxplot=TRUE, title="Sample boxplot from metabolite values") {
  
  require(reshape)
  
  # Generate new sample id
  data$Samples <- paste(id.col.abbrev, data[,id.colname],                 
                           group.a.abbrev, data[,group.a.colname],                   
                           group.b.abbrev, data[,group.b.colname],                   
                           sep="")
  data <- na.omit(data[c("Samples", id.colname, group.a.colname, group.b.colname, responses)])
  x.col <- "Samples"
  
  # Melt data for ggplot to make boxplots from multiple responses
  data <- melt(data,
               id.vars=c(x.col, group.a.colname, group.b.colname, id.colname),
               measure.vars=responses)
  # Order data by group.b.colname (usually TIME), group.a.colname (usually Gest Group), id.colname (usually subject/patient ID)
  data <- data[order(data[,group.b.colname], data[,group.a.colname], data[,id.colname]),]
  
  # Check levels or make them accordingly to designated column
  if(!is.null(x.levels)){
    data <- droplevels(data[data[,x.col] %in% x.levels,])
  } else { # is NULL
    # Make sure x.levels as factor and get levels for ordering
    x.levels <- levels(factor(data$Samples, levels=unique(as.character(data$Samples))))
  }
  
  data[,id.colname] <- as.factor(as.character(data[,id.colname]))
  
  # Plotting
  ###########
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  p <- ggplot(data, aes_string(x=x.col, y="value"))#, group=group.col, fill=group.col))
  # Colour boxplots by Group
  p <- p + geom_boxplot(aes_string(fill=group.a.colname))
  
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



draw.multiview.scatterplots <- function(data, file, responses, phenotypes, phenotypes.labels, x.levels=NULL, limit.amount=NULL, group.col=NULL, group.col.label=NULL, group.levels=NULL, group.stat_smooth=FALSE, facet.formula=NULL, colour="black", per.page.cols=length(phenotypes), per.page.rows=1, pdf.width=11, pdf.height=6) {
  
  require(MASS)
  require(splines)
  require(gridExtra)
  
  tmp.data <- data
  
  plots <- list()
  
  if(is.null(limit.amount)){
    limit.amount <- length(responses)
  }
  
  # Multicore, !order plots to response order will change
  ## create a list of plots, two different phenotype plots sequenzially
  plots <- foreach(i=1:limit.amount, .combine=c) %dopar% { 
    #for(i in 1:limit.amount){
    
    response <- responses[i]
    
    # Generate plot for each phenotype
    phenos.plots <- list()
    
    # %do% is sequential
    plots2 <- foreach(j=1:length(phenotypes), .combine=c) %do% {
      phenotype <- phenotypes[j]
      phenotype.label <- phenotypes.labels[j]
      
      # Select only data that matches group.levels in column group.col
      if(!is.null(group.levels)){
        tmp.data <- droplevels(tmp.data[tmp.data[,group.col] %in% group.levels,])
      } else if(!is.null(group.col)) { # is NULL
        group.levels <- levels(data[,group.col])
      }
      # set x.levels if not set
      if(!is.null(x.levels)){
        tmp.data <- droplevels(tmp.data[tmp.data[,phenotype] %in% x.levels,])
      } else { # is NULL
        x.levels <- levels(data[,phenotype])
      }
      
      p1 <- draw.scatterplot(data, response, phenotype, phenotype.label, x.levels, group.col, group.col.label, group.levels, group.stat_smooth,  colour, facet.formula)
      list(p1)
    }
    plots2
  }
  
  # Arrange plots, according to parametres or side by side (1 row, as many cols as phenotypes)
  p <- do.call(marrangeGrob, c(plots, ncol=per.page.cols, nrow=per.page.rows))
  
  tryCatch({
    ggsave(file, p, width=pdf.width, height=pdf.height)
    log.text(paste("Saved multiview scatter plots to: ", output.file, sep=""), log.file=log.file, logging=logging)
  }, error = function(e) {
    print(e$message)
    log.text(paste("Error, could not save the following file: ", output.file, sep=""), log.file=log.file, logging=logging)
  })
}

draw.scatterplot <- function(data, response, phenotype, phenotype.label, x.levels=NULL, group.col=NULL, group.col.label=NULL, group.levels=NULL, group.stat_smooth=FALSE,  colour="black", facet.formula=NULL) {
  
  tmp.data <- data 
  
  # Plotting
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  # Set to factor if only few unique data points
  if(length(unique(tmp.data[,phenotype])) < 5){
    tmp.data[,phenotype] <- factor(tmp.data[,phenotype])
    # Get limits
    x.limits <- levels(tmp.data[,phenotype])
  }
  
  p <- ggplot(tmp.data, aes_string(x=phenotype, y=response))#, group=group.col, fill=group.col))
  # p <- p + stat_smooth(colour=colour, method="loess")#, aes_string(group=group.col, colour=group.col))
  # Draw stat_smooth loess per group
  if(!is.null(group.col) & group.stat_smooth == TRUE){
    p <- p + stat_smooth(method=lm, family="binomial", formula=y~ns(x,2), aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
    p <- p + guides(fill=FALSE)
  }
  
  p <- p + geom_point(aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col), size=2) +
    labs(
      # Legend title for group
      group=group.col.label, fill=group.col.label, shape=group.col.label, colour=group.col.label)
  
  # Set discrete (With limits) or continuous x axis (with stat_smooth)
  if(length(unique(tmp.data[,phenotype])) < 5){
    p <- p + scale_x_discrete(limits=x.limits)
  } else {
    p <- p + scale_x_continuous(phenotype.label)
    if(is.null(group.col) & group.stat_smooth==FALSE){
      p <- p + stat_smooth(colour=colour, method="loess", aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
    }
  }
  
  p <- p + theme(legend.position="top") + xlab(phenotype.label)
  p <- p + ggtitle(response) + 
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, hjust=0.5))
  p <- p + theme_bw()
  
  p
}


draw.multiview.scatterplots <- function(data, file, responses, phenotypes, phenotypes.labels, x.levels=NULL, limit.amount=NULL, group.col=NULL, group.col.label=NULL, group.levels=NULL, group.stat_smooth=FALSE, colour="black", per.page.cols=length(phenotypes), per.page.rows=1, pdf.width=11, pdf.height=6) {
  
  require(MASS)
  require(splines)
  require(gridExtra)
  
  tmp.data <- data
  
  plots <- list()
  
  if(is.null(limit.amount)){
    limit.amount <- length(responses)
  }
  
  # Multicore, !order plots to response order will change
  ## create a list of plots, two different phenotype plots sequenzially
  plots <- foreach(i=1:limit.amount, .combine=c) %dopar% { 
    #for(i in 1:limit.amount){
    
    response <- responses[i]
    
    # Generate plot for each phenotype
    phenos.plots <- list()
    
    # %do% is sequential
    plots2 <- foreach(j=1:length(phenotypes), .combine=c) %do% {
      phenotype <- phenotypes[j]
      phenotype.label <- phenotypes.labels[j]
      
      # Select only data that matches group.levels in column group.col
      if(!is.null(group.levels)){
        tmp.data <- droplevels(tmp.data[tmp.data[,group.col] %in% group.levels,])
      } else if(!is.null(group.col)) { # is NULL
        group.levels <- levels(data[,group.col])
      }
      # set x.levels if not set
      if(!is.null(x.levels)){
        tmp.data <- droplevels(tmp.data[tmp.data[,phenotype] %in% x.levels,])
      } else { # is NULL
        x.levels <- levels(data[,phenotype])
      }
      
      p1 <- draw.scatterplot(data, response, phenotype, phenotype.label, x.levels, group.col, group.col.label, group.levels, group.stat_smooth,  colour)
      list(p1)
    }
    plots2
  }
  
  # Arrange plots, according to parametres or side by side (1 row, as many cols as phenotypes)
  p <- do.call(marrangeGrob, c(plots, ncol=per.page.cols, nrow=per.page.rows))
  
  tryCatch({
    ggsave(file, p, width=pdf.width, height=pdf.height)
    log.text(paste("Saved multiview scatter plots to: ", output.path, sep=""), log.file=log.file, logging=logging)
  }, error = function(e) {
    print(e$message)
    log.text(paste("Error, could not save the following file: ", output.file, sep=""), log.file=log.file, logging=logging)
  })
}

draw.scatterplot <- function(data, response, phenotype, phenotype.label, x.levels=NULL, group.col=NULL, group.col.label=NULL, group.levels=NULL, group.stat_smooth=FALSE,  colour="black") {
  
  tmp.data <- data 
  
  # Plotting
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  # Set to factor if only few unique data points
  if(length(unique(tmp.data[,phenotype])) < 5){
    tmp.data[,phenotype] <- factor(tmp.data[,phenotype])
    # Get limits
    x.limits <- levels(tmp.data[,phenotype])
  }
  
  p <- ggplot(tmp.data, aes_string(x=phenotype, y=response))#, group=group.col, fill=group.col))
  # p <- p + stat_smooth(colour=colour, method="loess")#, aes_string(group=group.col, colour=group.col))
  # Draw stat_smooth loess per group
  if(!is.null(group.col) & group.stat_smooth == TRUE){
    p <- p + stat_smooth(method=lm, family="binomial", formula=y~ns(x,2), aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
    p <- p + guides(fill=FALSE)
  }
  
  p <- p + geom_point(aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col), size=2) +
    labs(
      # Legend title for group
      group=group.col.label, fill=group.col.label, shape=group.col.label, colour=group.col.label)
  
  # Set discrete (With limits) or continuous x axis (with stat_smooth)
  if(length(unique(tmp.data[,phenotype])) < 5){
    p <- p + scale_x_discrete(limits=x.limits)
  } else {
    p <- p + scale_x_continuous(phenotype.label)
    if(is.null(group.col) & group.stat_smooth==FALSE){
      p <- p + stat_smooth(colour=colour, method="loess", aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
    }
  }
  
  p <- p + theme(legend.position="top") + xlab(phenotype.label)
  p <- p + ggtitle(response) + 
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, hjust=0.5))
  p <- p + theme_bw()
  
  p
}


# 
#   # Arrange plots, according to parametres or side by side (1 row, as many cols as phenotypes)
#   p <- do.call(marrangeGrob, c(plots, ncol=per.page.cols, nrow=per.page.rows))
#   
#   tryCatch({
#     ggsave(file, p, width=pdf.width, height=pdf.height)
#     log.text(paste("Saved multiview scatter plots to: ", output.path, sep=""), log.file=log.file, logging=logging)
#   }, error = function(e) {
#     print(e$message)
#     log.text(paste("Error, could not save the following file: ", output.file, sep=""), log.file=log.file, logging=logging)
#   })
# }

draw.scatterplot <- function(data, response, phenotype, phenotype.label, x.levels=NULL, group.col=NULL, group.col.label=NULL, group.levels=NULL, group.stat_smooth=FALSE,  colour="black") {
  
  tmp.data <- data 
  
  # Plotting
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  # Set to factor if only few unique data points
  if(length(unique(tmp.data[,phenotype])) < 5){
    tmp.data[,phenotype] <- factor(tmp.data[,phenotype])
    # Get limits
    x.limits <- levels(tmp.data[,phenotype])
  }
  
  p <- ggplot(tmp.data, aes_string(x=phenotype, y=response))#, group=group.col, fill=group.col))
  # p <- p + stat_smooth(colour=colour, method="loess")#, aes_string(group=group.col, colour=group.col))
  # Draw stat_smooth loess per group
  if(!is.null(group.col) & group.stat_smooth == TRUE){
    p <- p + stat_smooth(method=lm, family="binomial", formula=y~ns(x,2), aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
    p <- p + guides(fill=FALSE)
  }
  
  p <- p + geom_point(aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col), size=2) +
    labs(
      # Legend title for group
      group=group.col.label, fill=group.col.label, shape=group.col.label, colour=group.col.label)
  
  # Set discrete (With limits) or continuous x axis (with stat_smooth)
  if(length(unique(tmp.data[,phenotype])) < 5){
    p <- p + scale_x_discrete(limits=x.limits)
  } else {
    p <- p + scale_x_continuous(phenotype.label)
    if(is.null(group.col) & group.stat_smooth==FALSE){
      p <- p + stat_smooth(colour=colour, method="loess", aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
    }
  }
  
  # Theme, legend, xlabel, title
  p <- p + theme(legend.position="top") + xlab(phenotype.label)
  p <- p + ggtitle(response) + 
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, hjust=0.5))
  p <- p + theme_bw()
  
  p
}

draw.scatterplots <- function(data, responses, phenotypes, phenotypes.labels, x.levels=NULL, limit.amount=NULL, group.col=NULL, group.col.label=NULL, group.levels=NULL, group.stat_smooth=FALSE,  colour="black") {
  
  require(MASS)
  require(splines)
  tmp.data <- data
  
  phenos.plots <- list()
  for(i in 1:length(phenotypes)){
    
    phenotype <- phenotypes[i]
    phenotype.label <- phenotypes.labels[i]
    
    # Select only data that matches group.levels in column group.col
    if(!is.null(group.levels)){
      tmp.data <- droplevels(tmp.data[tmp.data[,group.col] %in% group.levels,])
    } else if(!is.null(group.col)) { # is NULL
      group.levels <- levels(data[,group.col])
    }
    
    if(!is.null(x.levels)){
      tmp.data <- droplevels(tmp.data[tmp.data[,phenotype] %in% x.levels,])
    } else { # is NULL
      x.levels <- levels(data[,phenotype])
    }
    
    # Plotting
    cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
    
    plots <- list()
    
    if(is.null(limit.amount)){
      limit.amount <- length(responses)
    }
    
    plots <- foreach(i=1:limit.amount, .combine=c) %dopar% { 
      library(MASS)
      library(splines)
      library(ggplot2)
      
      response <- responses[i]
      # Set to factor if only few unique data points
      if(length(unique(tmp.data[,phenotype])) < 5){
        tmp.data[,phenotype] <- factor(tmp.data[,phenotype])
        # Get limits
        x.limits <- levels(tmp.data[,phenotype])
      }
      
      p <- ggplot(tmp.data, aes_string(x=phenotype, y=response))#, group=group.col, fill=group.col))

      # p <- p + stat_smooth(colour=colour, method="loess")#, aes_string(group=group.col, colour=group.col))
      # Draw stat_smooth loess per group
      if(!is.null(group.col) & group.stat_smooth == TRUE){
        p <- p + stat_smooth(method=lm, family="binomial", formula=y~ns(x,2), aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
        p <- p + guides(fill=FALSE)
      }
      
      p <- p + geom_point(aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col), size=2) +
        labs(
          # Legend title for group
          group=group.col.label, fill=group.col.label, shape=group.col.label, colour=group.col.label)
      #+
       # guides(fill=guide_legend(title=group.col.label))
      
      # Draw line between ID /individual 
#       if(!is.null()){
#         p <- p + geom_line(aes_string(group="SUBJECT_ID"), alpha=0.1)
#       }
      
      # Set discrete (With limits) or continuous x axis (with stat_smooth)
      if(length(unique(tmp.data[,phenotype])) < 5){
        p <- p + scale_x_discrete(limits=x.limits)
      } else {
        p <- p + scale_x_continuous(phenotype.label)
        if(is.null(group.col) & group.stat_smooth==FALSE){
          p <- p + stat_smooth(colour=colour, method="loess", aes_string(group=group.col, fill=group.col, shape=group.col, colour=group.col))
        }
      }

      p <- p + theme(legend.position="top") + xlab(phenotype.label)
      p <- p + ggtitle(response) + 
        theme(plot.title = element_text(face="bold"),
              axis.text.x = element_text(angle=90, hjust=0.5))
      p <- p + theme_bw()
      
      list(p)
    }
    phenos.plots <- append(phenos.plots, plots)
  }
  phenos.plots
}

# Extract order of metabolites from metabolite subset by hierarchical clustering,
# create a combined group variable from Gest_group and Timepoint,
# aggregate data by function (default is mean,)
# draw heatmap from aggregated data
# Usage example:
#   mg_heatmap <- metabolite.grouping.heatmap(data=data, responses=responses,
#       group.a.colname="PREBABY_BACKGROUND_Gest_Group", group.a.abbrev="G"
#       group.b.colname="PREBABY_BACKGROUND_Time_code", group.b.abbrev="T"
#       combined.group.label="Group_and_Time")
metabolite.grouping.heatmap <- function(data, responses,
                                        group.a.colname, group.a.abbrev,
                                        group.b.colname, group.b.abbrev,
                                        combined.group.label,
                                        fun="mean",
                                        colour.low="white",
                                        colour.mid=NULL,
                                        colour.high="steelblue",
                                        x.text.angle=0){
  require(reshape)
  
  tmp.data <- data
  
  # create new variable with group and time combined
  tmp.data$combined <- paste(group.a.abbrev, data[,group.a.colname,],
                             group.b.abbrev, data[,group.b.colname,], sep="")
  # make it a factor
  tmp.data$combined <- as.factor(tmp.data$combined)
  tmp.data <- na.omit(tmp.data[c("combined", group.a.colname, group.b.colname, responses)])
  tmp.data <- droplevels(tmp.data)
  
  data.dist <- dist(t(tmp.data[responses]))
  #data.dist <- dist(data[responses])
  hc <- hclust(data.dist)
  hc.order <- hc$labels[c(hc$order)]
  
  # Aggregate data by mean of the metabolite per combined group and time factor (GT), with function mean as default
  melted <- melt.responses.by.groups(tmp.data, responses, "combined", x.label=combined.group.label, y.label="Metabolite", fun=fun)
  # reorder Metabolites to hclust order
  melted$Metabolite <- factor(melted$Metabolite, levels=hc.order, ordered=TRUE)  
  
  p1 <- draw.heatmap(melted, x.label=combined.group.label, y.label="Metabolite", x.axis.text.angle=x.text.angle, colour.low=colour.low, colour.mid=colour.mid, colour.high=colour.high)
  p1
  
}


# Create new variable from subject, group and time as a sample,
# count distances between samples (from all the subset metabolites),
# extract sample order by hierarchical clustering,
# draw sample vs sample heatmap from distance between them (where sample order is from hiearchical clustering)

# <- samples.heatmap(data=data, responses=responses,
#           id.colname="SUBJECT_ID", id.col.abbrev="SUBJ",
#           group.a.colname="meal", group.a.abbrev="_G",
#           group.b.colname="time", group.b,abbrev="_T")
samples.heatmap <- function(data, responses,
                            id.colname, id.col.abbrev,
                            group.a.colname, group.a.abbrev,
                            group.b.colname, group.b.abbrev,
                            separate.group.bar, separate.group.bar.label,
                            dist.method="euclidean",
                            colour.low="white",
                            colour.mid=NULL,
                            colour.high="steelblue"){
  require(reshape)
  require(gridExtra)
  
  tmp.data <- data
  
  # Generate new sample id
  tmp.data$new.id <- paste(id.col.abbrev, data[,id.colname],
                           group.a.abbrev, data[,group.a.colname],
                           group.b.abbrev, data[,group.b.colname],
                           sep="")
  
  tmp.data <- na.omit(tmp.data[c("new.id", id.colname, group.a.colname, group.b.colname, responses)])
  # Sample "names"
  #print(tmp.data$new.id)
  rownames(tmp.data) <- tmp.data$new.id
  
  sample.distances <- dist(tmp.data[responses], method=dist.method)
  # Get order by hiearchical clustering 
  hc <- hclust(sample.distances)
  hc.order <- hc$labels[c(hc$order)]
  
  distances <- as.matrix(sample.distances)
  melted <- melt(distances)
  melted<- merge(melted, tmp.data[c("new.id", separate.group.bar)], by.x="X1", by.y="new.id")
  
#   testi <- as.data.frame(testi)
#   testi <- cbind(testi, group=tmp.data$Group)
#   
#   melted <- melt(testi, id.vars=c("group"))
#   
#   testi.2 <- melt(testi)
#   
#   testi.2$new.id <- factor(testi.2$new.id, levels=hc.order, ordered=TRUE)
#   testi.2$X2 <- factor(testi.2$X2, levels=hc.order, ordered=TRUE)
#   
#   testi.2 <- rename(testi.2, (c("new.id"="X1")))
#   melted <- melt(testi, id.vars=c("group"), variab)

  # reorder Samples to hclust order
  melted$X1 <- factor(melted$X1, levels=hc.order, ordered=TRUE)
  melted$X2 <- factor(melted$X2, levels=hc.order, ordered=TRUE) 
  
  p <- draw.heatmap(melted, x.axis.text.angle=90,
                    title=expression(atop("Heatmap", atop(italic("Samples ordered by hierarchical clustering. Euclidean distance between samples."), ""))),
                    colour.low=colour.low,
                    colour.mid=colour.mid,
                    colour.high=colour.high)

  if(!is.null(separate.group.bar)){  
    # draw coloured bar
    p0 <- ggplot(data=melted, aes_string(x="X1", y=1)) +
      geom_tile(aes_string(fill=separate.group.bar)) + 
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks=element_blank(),
            axis.text=element_blank())

    # scale_y_discrete(limits=0.01)
    # scale_y_continuous(expand=c(0,0))
  
    gp1 <- ggplotGrob(p)
    gp2 <- ggplotGrob(p0)
    
    maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
    gp1$widths[2:5] <- as.list(maxWidth)
    gp2$widths[2:5] <- as.list(maxWidth)
    
    
    p <- arrangeGrob(gp1,gp2, ncol=1, heights=c(10/11,1/11))
    
#     pdf(plot.file, width=18, height=20)
#     grid.arrange(p)
#     dev.off()
#     
#     ########################
#     plot.file <- paste(output.prefix, "heatmap_samples_vs_samples.pdf", sep="")
#     ggsave(file=plot.file, scale=7, limitsize=FALSE)
#     grid.arrange(p, ncol=1, heights=c(10/11,1/11))
#     dev.off()
#     
#     p <- arrangeGrob(gp1,gp2, ncol=1)
#     
#     ########################
#     pdf(plot.file, pointsize=20)
#     p <- grid.arrange(gp1,gp2, ncol=1, heights=c(10/11,1/11))
#     dev.off()
  }
  
  p
}


# 1. Calculate distance between data variables
# 2. Do hieararchical clustering
# 3. Plot a dendrogram
hclust.dendro <- function(data, responses, transpose.data=FALSE, title="Hierarchical clustering", dist.method="euclidean", hclust.method="complete", rotate=TRUE, size=4, theme_dendro=FALSE) {
  require(ggplot2)
  require(ggdendro)
  
  data <- data[responses]
  
  if(transpose.data){
    data <- t(data)
  }
  
  dist.data <- dist(data, method=dist.method)
  hclust.data <- hclust(dist.data, method=hclust.method)
  
  p <- ggdendrogram(
    hclust.data, 
    rotate = rotate, 
    size = size, 
    theme_dendro = theme_dendro
  ) + labs(title=title)
  
  p
  
  # FOR GGPLOT
  #   hcdata <- dendro_data(hclust.data, type="rectangle")
  #   ggplot() + 
  #     geom_segment(data=segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
  #     geom_text(data=label(hcdata), aes(x=x, y=y, label=label, hjust=0), size=3) +
  #     coord_flip() + scale_y_reverse(expand=c(0.2, 0))
  
}


select.prefixed.cols <- function(data, prefix.list) {
  x <- data.frame()
  for (prefix in prefix.list){
    if(length(x) > 0){
      x = cbind(x, data[,grep(prefix, colnames(data))])
    }
    else{ 
      x = data[,grep(prefix, colnames(data))] 
    }
  }
  x
}



## Aggregate data by selected group.levels from a selection of variables(responses) with function (mean by default)
## Then rescale values from 0 to 1.
melt.responses.by.groups <- function(data, responses, group.col, group.levels=NULL, x.col=NULL, x.levels=NULL, fun="mean", x.label, y.label){
   
  # for melt
  require(reshape)
  require(scales)
  require(plyr)
  
  ## Adjusting data
  tmp.data <- data
  
  # Select data according to group.levels
  if(!is.null(group.levels)){
    tmp.data <- droplevels(tmp.data[tmp.data[,group.col] %in% group.levels,])
  } else { # is NULL
    group.levels <- levels(as.factor(tmp.data[,group.col]))
  }
  if(!is.null(x.col)){
    if(!is.null(x.levels)){
      tmp.data <- droplevels(tmp.data[tmp.data[,x.col] %in% x.levels,])
    } else { # is NULL
      x.levels <- levels(as.factor(tmp.data[,x.col]))
    }
  }
  
  # Check and log if data contains NA's
  if(any(is.na(tmp.data[responses]))){
    # "Data contains NA values! 'na.rm' is set true"
  }
  
  # Select wanted responses according to groups.levels and compute by function
  # Multigroup for example time and meal group
  if(!is.null(x.col)){
    # If we want to leave also another column defined by x.col (time for example)
    tmp.data <- aggregate(tmp.data[responses], by=list(a=tmp.data[, group.col], b=tmp.data[, x.col]), FUN=fun, na.rm=TRUE)
    tmp.data <- melt(tmp.data, id.vars=c("a", "b")) 
    tmp.data <- ddply(tmp.data, .(variable), transform, rescale=rescale(value))
    # rename columns
    tmp.data <- rename(tmp.data, (c("a"=x.label, "b"=x.col, "variable"=y.label)))
  }
  else{
    tmp.data <- aggregate(tmp.data[responses], by=list(a=tmp.data[, group.col]), FUN=fun, na.rm=TRUE)
    tmp.data <- melt(tmp.data, id.var="a")
    tmp.data <- ddply(tmp.data, .(variable), transform, rescale=rescale(value))
    # rename columns
    tmp.data <- rename(tmp.data, (c("a"=x.label, "variable"=y.label)))
  }
  tmp.data
}


# Compute correlation
compute.corr <- function(data, var1=NULL, var2=NULL, use="everything", method="spearman"){
  
  if(is.null(var1)){
    # Correlation all vs all
    cor.data <- cor(data, use=use, method=method)
  }
  else{
    if(is.null(var2))
      ## TODO: var1 variables vs each other or all???
      cor.data <- cor(data[,var1], use=use, method=method)
    else{
      cor.data <- cor(data[,var1], data[,var2], use=use, method=method)
    }
  }
  
## TODO: return melted???? (better for draw.heatmap)
  # melt(cor.data)
  cor.data
}


# NEEDS REVIEW/CHECK!
plot.tsne <- function(data, grouping, grouping.title, perp=50, max.iter=1000, title="t-SNE" ){
  
  require(tsne)
  
  # tsne example
#   colors = rainbow(length(unique(iris$Species)))
#   names(colors) = unique(iris$Species)
#   ecb = function(x,y){ plot(x,t='n'); text(x,labels=iris$Species, col=colors[iris$Species]) }
#   tsne_fef = tsne(iris[,1:4], epoch_callback = ecb, perplexity=50)

  # make variables accessible outside/ inside ecb function
  grouping <<- grouping
  grouping.title <<- grouping.title

  # Graphing function for tsne
  ecb = function(x){
    
    # because ggplot
    x = data.frame(x)
    names(x) = c('x', 'y')
    
    p <- ggplot(x, aes(x,y))
    p <- p + geom_point(aes(colour=grouping,shape=grouping),size=2.5) +
      labs(
       # Legend title for group
       color=grouping.title, shape=grouping.title)
    p <- p + ggtitle(title)
  
    #print(p)
  }
  
  # run tsne
  tsne_ = tsne(
    data,
    epoch_callback = ecb,
    perplexity=perp,
    max_iter=max.iter
  )

  tsne_
}


compute.pca <- function(data, pca.scale=FALSE, pca.center= FALSE){
  
  pca <- prcomp(data, scale.=pca.scale, center=pca.center)
 
  pca
}

plot.pca.gradient <- function(pca.data, first.group, first.group.title, component.x = 1, component.y = 2, title="PCA", sampleid.labels = NULL, theme=NULL, facet.group=NULL, color.col="ids", gradient=FALSE, text.labels=TRUE, fixed_coordinates=FALSE){
  
  require(MASS)
  require(ggplot2)
  require(scales)
  #require(gridExtra)
  
  # Capture local enviroment, because problematic ggplot aes with local variables
  .e <- environment()
  
  # Grouping variable and pca results to same dataframe
  scores = data.frame(group=first.group, pca=pca.data$x)
  scores$group = factor(scores$group)
  if(!is.null(sampleid.labels)){
    scores <- cbind(scores, ids=sampleid.labels)
  }
  if(!is.null(facet.group)){
    scores <- cbind(scores, facet=as.factor(facet.group))
  }
  
  # for percentage labels
  prop.pca = pca.data$sdev^2/sum(pca.data$sdev^2)
  
  # string variables made from component numbers to be referenced 
  PCX <- paste("pca.PC", as.character(component.x), sep="")
  PCY <- paste("pca.PC", as.character(component.y), sep="")
  
  p <- ggplot(scores, environment=.e) +
    geom_point(aes_string(x=PCX, y=PCY, colour=color.col, shape="group"),size=2.5) +
    labs(x = paste(PCX, "(", percent(prop.pca[component.x]), ")", sep=""),
         y = paste(PCY, "(", percent(prop.pca[component.y]), ")", sep=""),
         # Legend title for group
         shape=first.group.title)
  
  # Show datapoint(sample/id) labels if not null
  if(!is.null(sampleid.labels)){
    if(text.labels){
      p <- p + geom_text(aes_string(x=PCX, y=PCY, label="ids"), size=4, hjust=1.5)
    }
    
    if(gradient){
      p <- p + scale_colour_gradient(low="yellow",high="red")
    }
  }
  
  
  if(!is.null(facet.group)){
    
    x <- paste(". ~", "facet")
    p <- p + facet_grid(. ~ facet)#, labeller=label_both) #p <- p + opts(legend.position = "none")
    
  }
  #p <- p + geom_line(aes_string(x=PCX, y=PCY, group="ids", colour="facet"))
  
  # X,Y-axis to same aspect ratio
  if(fixed_coordinates){ p <- p + coord_fixed() }
  
  # Apply theme
  if(!is.null(theme)){ # Fill 
  }
  else{ p <- p + theme_bw() }
  
  #p <- p + theme(legend.position="top")
  p <- p + ggtitle(title)
  p
}

plot.pca <- function(pca.data, first.group, first.group.title, component.x = 1, component.y = 2, title="PCA", sampleid.labels = NULL, theme=NULL, facet.group=NULL, fixed_coordinates=FALSE){
  
  require(MASS)
  require(ggplot2)
  require(scales)
  require(ggrepel)
  #require(gridExtra)
  
  # Capture local enviroment, because problematic ggplot aes with local variables
  .e <- environment()
  
  # Grouping variable and pca results to same dataframe
  if(!is.null(first.group)){
    scores = data.frame(group=first.group, pca=pca.data$x)
    scores$group = factor(scores$group)
  } else{
    scores = data.frame(pca=pca.data$x)
  }
  if(!is.null(sampleid.labels)){
    scores <- cbind(scores, ids=sampleid.labels)
  }
  if(!is.null(facet.group)){
    scores <- cbind(scores, facet=as.factor(facet.group))
  }
  
  # for percentage labels
  prop.pca = pca.data$sdev^2/sum(pca.data$sdev^2)
  
  # string variables made from component numbers to be referenced 
  PCX <- paste("pca.PC", as.character(component.x), sep="")
  PCY <- paste("pca.PC", as.character(component.y), sep="")
  
  if(is.null(first.group.title)){
    p <- ggplot(scores, environment=.e) +
      geom_point(aes_string(x=PCX, y=PCY),size=2.5) +
      labs(x = paste(PCX, "(", percent(prop.pca[component.x]), ")", sep=""),
           y = paste(PCY, "(", percent(prop.pca[component.y]), ")", sep=""))
          # Legend title for group
  }
  
  # Color and shape by group if exists
  if(!is.null(first.group.title)){
    p <- ggplot(scores, environment=.e) +
      geom_point(aes_string(x=PCX, y=PCY, colour="group", shape="group"),size=2.5) +
      labs(x = paste(PCX, "(", percent(prop.pca[component.x]), ")", sep=""),
           y = paste(PCY, "(", percent(prop.pca[component.y]), ")", sep=""),
           color=first.group.title, shape=first.group.title)
  }
  
  # Show datapoint(sample/id) labels if not null
  if(!is.null(sampleid.labels)){
    p <- p + geom_text_repel(aes_string(x=PCX, y=PCY, label="ids"), size=4)
  }
  
  # X,Y-axis to same aspect ratio
  if(fixed_coordinates){ p <- p + coord_fixed() }
  
  if(!is.null(facet.group)){
    
    x <- paste(". ~", "facet")
    p <- p + facet_grid(. ~ facet)#, labeller=label_both) #p <- p + opts(legend.position = "none")
    
  }
  #p <- p + geom_line(aes_string(x=PCX, y=PCY, group="ids", colour="facet"))
  
  
  # Apply theme
  if(!is.null(theme)){ # Fill 
  }
  else{ p <- p + theme_bw() }
  
  #p <- p + theme(legend.position="top")
  p <- p + ggtitle(title)
  p
}

# Subset data according to group.a,
# calculate PCA
# Colour results according to group.b
plot.multigroup.pca <- function(data, responses, group.subset, group.subset.title, group.coloring, group.coloring.title, title.prefix="PCA", sampleid.labels= NULL, grid.row=NULL, grid.col=NULL, scale=FALSE, center=FALSE, fixed_coordinates=FALSE){
  
  require(gridExtra)
  
  plots <- list()
  
  if(!is.factor(data[,group.subset])){
    data[,group.subset] <- factor(data[,group.subset])
  }
  if(is.null(grid.col) & is.null(grid.row)){
    grid.col <- length(levels(data[,group.subset]))
    grid.row <- 1
  }
  
  for(x in levels(data[,group.subset])){
    id.labels=NULL
    # Subset data by group.subset for PCA computing
    data.from.group <- data[data[,group.subset] == x,][responses]
    # Subset other group coloring to match previous subset 
    group.coloring.data <- data[data[,group.subset] == x,][group.coloring]
    # Subset sample/id labels to match previous subset
    if(!is.null(sampleid.labels)){
      id.labels <- data[data[,group.subset] == x,][sampleid.labels]
    }
    #print(sampleid.labels)
    # Compute PCA for subset
    pca.results <- compute.pca(data.from.group, pca.scale=scale, pca.center=center)
    # Rename title according subset group
    new.title <- paste(title.prefix, group.subset.title, x, sep=" ")
    
    # Warning and skipping plotting for pca.results when there is only fewer than 2 principal components generated
    if(ncol(pca.results$x) < 2){
      plots <- plots
      info.string <- paste("Warning, check your data! Value",x ,"in", group.subset, "ONLY produced PC1 in PCA. Skipping plot generation for", group.subset,"value", x, sep=" ")
      warning(info.string)
    } else {
      p <- plot.pca(pca.results, group.coloring.data[,group.coloring], group.coloring.title, title=new.title, sampleid.labels=id.labels[,sampleid.labels])
      
      plots <- c(plots, list(p))
      rm(p)
    }

  }
  #do.call("grid.arrange", c(plots, ncol=grid.col, nrow=grid.row))
  plots  
}


compute.pls_da <- function(a, b, n.comp=2){
  
  require("mixOmics")
  
  result <- plsda(a, b, n.comp)
  result  
}


plot.pls_da <- function(data, sample.groups, group.title, component.x=1, component.y=2, title="PLS-DA"){

  require("ggplot2")
  #qplot(data$variates$X[,1], data$variates$X[,2], col=sample.groups))
  
  data.df <- as.data.frame(data$variates$X)
  data.x <- data$variates$X[,component.x]
  data.y <- data$variates$X[,component.y]
  .e <- environment()
  
  p <- ggplot(data.df, aes(x = data.x, y = data.y, group=sample.groups, fill=sample.groups, colour=sample.groups), environment=.e)
  p <- p + geom_point() + labs(color=group.title, fill=group.title)
  
  p <- p + theme(legend.position="right")
  
  p <- p + ggtitle(title) #+ 
    
#     theme(plot.title = element_text(face="bold"),
#           axis.text.x = element_text(angle=90, hjust=0.5))
#   p <- p + theme_bw() 
#   
#   
#   p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)
#   #p <- p + scale_x_continuous(breaks=x.levels, labels=x.levels, limits=x.levels)
#   print(p)

  p
}


## Draw heatmap from aggregated data
## Draw tile fill by fill.column
## !!! Takes MELT type data
draw.heatmap <- function(data,
                         x.label="X1",
                         y.label="X2",
                         x.levels=NULL,
                         z.label=NULL, 
                         fill.column="rescale",
                         x.axis.text.angle=0,
                         value.limits=c(0,1),
                         md=NULL,
                         colour.low="white",
                         colour.mid=NULL,
                         colour.high="steelblue",
                         remove.axis.labels=FALSE,
                         title="Heatmap",
                         theme=NULL,
                         gradient.colours=NULL,
                         gradient.breaks=NULL,
                         gradient.limits=NULL) {

  # For lazyness sake
  # check if rescale column exists, otherwise replace with "value"
  # if "value" column then get value limits from min and max
  if(fill.column=="rescale"){
    if(!"rescale" %in% colnames(data))
      fill.column <- "value"
    
      max.value <- max(data[,fill.column], na.rm=TRUE)
      min.value <- min(data[,fill.column], na.rm=TRUE)
    
      value.limits=c(min.value,max.value)
  } else {
    if(!is.factor(data[, fill.column])){
      max.value <- max(data[,fill.column], na.rm=TRUE)
      min.value <- min(data[,fill.column], na.rm=TRUE)
      
      value.limits=c(min.value,max.value)
    } else {
      info.string <- "Values are factors, no value.limits defined!"
      warning(info.string)
      value.limits <- NULL
    }
  }

  #print("Value limits: ")
  #print(value.limits)

  # set new md if md is null
  if(is.null(md)){
    md <- max.value/2
  }
  
  # reverse row order / y.axis labels from top to down
  # default is from first at the bottom of y.axis
  data[,y.label] = with(data, factor(get(y.label), levels=rev(levels(get(y.label)))))
  
  # plotting
  p <- ggplot(data, aes_string(x=x.label, y=y.label)) + geom_tile(aes_string(fill = fill.column), colour = "white")
  
  # Apply theme
  if(!is.null(theme)){
    # Fill
    print("testi") 
  }else{ p <- p + theme_bw() }
  
  if(!is.null(z.label)){
    # Multigroup data, create facet grid
    x <- paste(".~", z.label)
    p <- p + facet_grid(x) #p <- p + opts(legend.position = "none")
  }
  
  p <- p + theme(axis.text.x = element_text(angle=x.axis.text.angle, hjust=1, vjust=1))
  
#   po.nopanel <- theme(panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank())
#   p <- p + scale_fill_brewer(palette = "RdYlGn",name="Correlation") # RColorBrewer package
#   p <- p + po.nopanel
  #palette <- colorRampPalette(gradient.colours)
  #p <- p + scale_fill_manual(values=gradient.colours, limits=value.limits, breaks=gradient.brekas)

#   print(gradient.limits[1])
#   print(value.limits[1])
#   print(gradient.limits[2])
#   print(value.limits[2])

  # If no mid color defined
  if(is.null(colour.mid)){
    if(!is.null(gradient.colours)){ # gradient colors defined
      # Checking input limits (min & max values for plotting)
      if(is.null(gradient.limits)){ # no gradient.limits input, give value.limits as default
        gradient.limits <- value.limits
      } else { # gradient.limits defined, check if that are outside
        if(gradient.limits[1] > value.limits[1]){
          warning.string <- "Defined (lower) gradient.limit smaller than minimum value limit in data! Data lost from image!"
          warning(warning.string)
        }
        if(gradient.limits[2] < value.limits[2]){
          warning.string <- "Defined (upper) gradient.limit smaller than maximum value limit in data! Data lost from image!"
          warning(warning.string)
        }
      }
      
      print("min & max values (correlation:")
      print(value.limits)
      
      p <- p + scale_fill_gradient(colours=rev(gradient.colours), limits=gradient.limits, breaks=gradient.breaks)
      
    } else {
      p <- p + scale_fill_gradient(low=colour.low, high=colour.high, limits=value.limits)
    }
  }else{ # mid colour defined
      p <- p + scale_fill_gradient2(low=colour.low, midpoint=md, mid=colour.mid, high=colour.high, limits=value.limits)
  }
  
  # Order x-axis according to x.levels if not null
  if(is.null(x.levels)){
    x.levels <- levels(as.factor(data[,x.label]))
  }
  p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)
  
  # Add title to plot
  p <- p + ggtitle(title)
  
  if(remove.axis.labels){
    p <- p + theme(
         axis.title.x=element_blank(),
         axis.title.y=element_blank()) + ylab("") + xlab("")
  }
  p
}


draw.plots.org <- function(data, responses, group.col, group.levels=NULL, x.col, x.levels=NULL, x.label, y.label, limit.amount=NULL) {
  
  plots <- list()
  count <- 0
  
  if(is.null(limit.amount)){
    limit.amount <- length(responses)
  }
  
  # Select only data that matches group.levels in column group.col
  if(!is.null(group.levels)){
    data <- droplevels(data[data[,group.col] %in% group.levels,]) 
  } else { # is NULL
    group.levels <- levels(data[,group.col])
  }
  
  if(!is.null(x.levels)){
    data <- droplevels(data[data[,x.col] %in% x.levels,])  
  } else { # is NULL
    x.levels <- levels(data[,x.col])
  }

  # Plotting
  # Mean point with errorbars, multiple groups with dodge
  pd <- position_dodge(1)
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  for (response in responses) {
    count <- count + 1
    # limit how many metabolites are returned for quick debugging etc
    if (count > limit.amount) break
    
    p <- ggplot(data, aes_string(x=x.col, y=response, group=group.col, fill=group.col))
    
    p <- p + geom_line(stat="summary",stat_params=list(fun.y="mean"), aes_string(colour=group.col, linetype=group.col), position=pd, size=1)
    p <- p + geom_point(stat="summary",stat_params=list(fun.y="mean"), aes_string(colour=group.col, shape=group.col), position=pd, size=3)
    p <- p + stat_summary(fun.data = "mean_se", aes_string(colour=group.col), geom="errorbar",position=pd, width=.5)
     
    p <- p + theme(legend.position="top")
    p <- p + ggtitle(response) + 
      theme(plot.title = element_text(face="bold"),     
            axis.text.x = element_text(angle=90, hjust=0.5))  
    p <- p + theme_bw()    
    p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)
    
    # Add to list  
    plots <- c(plots, list(p))   
    rm(p)
  }
  plots
}

draw.plots <- function(data, responses, group.col, group.col.title=NULL, group.levels=NULL, x.col, x.levels=NULL, limit.amount=NULL, position.dodge=0.3) {

  tmp.data <- data
  
  # Select only data that matches group.levels in column group.col
  if(!is.null(group.levels)){
    tmp.data <- droplevels(tmp.data[tmp.data[,group.col] %in% group.levels,])
  } else { # is NULL
    group.levels <- levels(data[,group.col])
  }
  
  if(!is.null(x.levels)){
    tmp.data <- droplevels(tmp.data[tmp.data[,x.col] %in% x.levels,])
  } else { # is NULL
    x.levels <- levels(data[,x.col])
  }
  
  # Plotting
  
  # Mean point with errorbars, multiple groups with dodge 
  pd <- position_dodge(position.dodge)
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  plots <- list()
    
  if(is.null(limit.amount)){
    limit.amount <- length(responses)
  }
  
  if(is.null(group.col.title)){
    group.col.title <- group.col
  }
  #gc()
  
  plots <- foreach(i=1:limit.amount, .combine=c) %dopar% { 
    require(ggplot2)
    
    response <- responses[i]
    #print(response)
    #print(i)
   
    p <- ggplot(tmp.data, aes_string(x=x.col, y=response, group=group.col, colour=group.col)) +
      # Errorbars with solid lines
      stat_summary(fun.data = "mean_se",
                   geom="errorbar",position=pd, width=.5,
                   fun.ymin = function(x) mean(x) - sd(x), 
                   fun.ymax = function(x) mean(x) + sd(x)) +
      # Plot point to mean
      stat_summary(fun.y = mean,
                   geom = "point",
                   position=pd, size=3,
                  fun.ymin = min,
                  fun.ymax = function(x) mean(x) + sd(x),
                  aes_string(shape=group.col)) +
      # Line from mean to mean between for example timepoints
      stat_summary(fun.y = mean,
                   geom = "line",
                   position=pd, size=0.5,
                   fun.ymin = min,
                   fun.ymax = function(x) mean(x) + sd(x),
                   aes_string(linetype=group.col))

    p <- p + theme(legend.position="top")
    p <- p + ggtitle(response) + 
      theme(plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle=90, hjust=0.5))
    p <- p + theme_bw() 

    # Check if discrete
    if(is.factor(x.levels) & !(is.numeric(x.levels))){
      p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)
    }

    # Check if completely numeric then make continuous
    if(is.null(x.levels) & is.numeric(data[,x.col])){
      p <- p + scale_x_continuous(breaks=unique(data[,x.col]), labels=unique(data[,x.col]))
    }
    if(is.factor(x.levels) & !(is.numeric(x.levels))){
      p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)
    }

    p <- p + scale_colour_discrete(name=group.col.title) +
      scale_shape_discrete(name=group.col.title) +
      scale_linetype_discrete(name=group.col.title) +
      scale_fill_discrete(name=group.col.title)
    
    list(p)
  }

  plots
}

draw.plots.tofile <- function(data, responses, group.col, group.col.title=NULL, group.levels=NULL, x.col, x.levels=NULL, limit.amount=NULL, position.dodge=0.3, output.file, log.description="Saved line plots to: ", log.file, logging) {
  
  tmp.data <- data
  
  # Select only data that matches group.levels in column group.col
  if(!is.null(group.levels)){
    tmp.data <- droplevels(tmp.data[tmp.data[,group.col] %in% group.levels,])
  } else { # is NULL
    group.levels <- levels(data[,group.col])
  }
  
  if(!is.null(x.levels)){
    tmp.data <- droplevels(tmp.data[tmp.data[,x.col] %in% x.levels,])
  } else { # is NULL
    x.levels <- levels(data[,x.col])
  }
  
  # Plotting
  
  # Mean point with errorbars, multiple groups with dodge 
  pd <- position_dodge(position.dodge)
  cbgColourPalette <- scale_colour_manual(values=c("#0072B2", "#D55E00", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7"))
  
  plots <- list()
  
  if(is.null(limit.amount)){
    limit.amount <- length(responses)
  }
  
  if(is.null(group.col.title)){
    group.col.title <- group.col
  }
  
  plots <- foreach(i=1:limit.amount) %dopar% { 
    
    response <- responses[i]
    #print(response)
    #print(i)
    
    p <- ggplot(tmp.data, aes_string(x=x.col, y=response, group=group.col, colour=group.col)) +
      # Errorbars with solid lines
      stat_summary(fun.data = "mean_se",
                   geom="errorbar",position=pd, width=.5,
                   fun.ymin = function(x) mean(x) - sd(x), 
                   fun.ymax = function(x) mean(x) + sd(x)) +
      # Plot point to mean
      stat_summary(fun.y = mean,
                   geom = "point",
                   position=pd, size=3,
                   fun.ymin = min,
                   fun.ymax = function(x) mean(x) + sd(x),
                   aes_string(shape=group.col)) +
      # Line from mean to mean between for example timepoints
      stat_summary(fun.y = mean,
                   geom = "line",
                   position=pd, size=0.5,
                   fun.ymin = min,
                   fun.ymax = function(x) mean(x) + sd(x),
                   aes_string(linetype=group.col))
    
    p <- p + theme(legend.position="top")
    p <- p + ggtitle(response) + 
      theme(plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle=90, hjust=0.5))
    p <- p + theme_bw() 
    
    # Check if discrete
    if(is.factor(x.levels) & !(is.numeric(x.levels))){
      p <- p + scale_x_discrete(breaks=x.levels, labels=x.levels, limits=x.levels)
    }
    
    # Check if completely numeric then make continuous
    if(is.null(x.levels) & is.numeric(data[,x.col])){
      p <- p + scale_x_continuous(breaks=unique(data[,x.col]), labels=unique(data[,x.col]))
    }
    
    p <- p + scale_colour_discrete(name=group.col.title) +
      scale_shape_discrete(name=group.col.title) +
      scale_linetype_discrete(name=group.col.title) +
      scale_fill_discrete(name=group.col.title)
    
    # Save each page as a separate file (filename ends with number of the response)
    pdf(paste0(output.file, i))
    plot(p)
    dev.off()
    
  }
  
  # Stitch individual response pdf:s together into a single file
  files <- dir(pattern=output.file)
#  system2(command = command = "pdftk",args = c(shQuote(files), "cat output", shQuote(output.file)))

  
  log.text(paste(log.description, output.path, sep=""), log.file=log.file, logging=logging)
  
}



save.plotlist.pdf <- function(output.path, plotlist, log.description, log.file, logging) {
  require(gridExtra)
  require(doMC)
  
#   listlength <- length(plotlist)
#   print(listlength)
#   page.limit <- 1000
  
#   if(listlength > page.limit){
#     info.string <- paste0("Trying to save more than 1000 pages! Splitting pdf to multiple files, 1000 pages per pdf...")
#     warning(info.string)
#     log.text(paste0(info.string), log.file=log.file, logging=logging)
#     
#     # Paging per document .pdf
#     tmp.df <- data.frame(n=seq(1, ceiling(listlength/page.limit), 1))
#     tmp.df$start <- c(1, seq(page.limit+1, listlength, page.limit))
#     if(listlength == tail(seq(page.limit, listlength, page.limit), n=1)){
#       tmp.df$end <- seq(page.limit, listlength, page.limit)
#     } else {
#       tmp.df$end <- c(seq(page.limit, listlength, page.limit), listlength)  
#     }
#     # print(tmp.df)
#     tmp.path <- output.path
#     
#     for(filu in 1:max(tmp.df$n)){
#       # Add file number end of file name
#       tmp.string <- paste0("_", filu,"of", max(tmp.df$n), ".pdf")
#       tmp.path <- gsub(".pdf", tmp.string, output.path)
# #       print(tmp.df[filu, "start"])
# #       print(tmp.df[filu,"end"])
#       pdf(tmp.path, onefile=TRUE)
#       #ggsave(output.path)
#       # Save pages
#       for (x in tmp.df[filu, "start"]:tmp.df[filu,"end"]) {
#         is.numeric(x)
#         #print(x)
#         tryCatch({
#           plot(plotlist[[x]])
#         }, error = function(err){
#           print(paste("ERROR during plot(plotlis[[x]]): ", err))
#           print(paste("Problem with x: ", x))
#           stop()
#         })
#       }
#       dev.off()
#     }
#   } else {
  pdf(output.path, onefile=TRUE)
  #ggsave(output.path)
  for (x in 1:length(plotlist)) {
    plot(plotlist[[x]])
  }
  dev.off()
# }
  log.text(paste(log.description, output.path, sep=""), log.file=log.file, logging=logging)
}

ggsave.plot <- function(plot, plot.file, log.description, log.file, logging, scale=1, width=par("din")[1], height=par("din")[2]){
  require(ggplot2)
  ggsave(plot, file=plot.file, scale=scale, width=width, height=height, limitsize=FALSE)
  dev.off()
  log.text(paste(log.description, plot.file, sep=""), log.file=log.file, logging=logging)

}

save.plot <- function(plot, plot.file, log.description, log.file, logging, scale=1, width=par("din")[1], height=par("din")[2]){
  pdf(plot, file=plot.file, width=width, height=height)
  plot(plot)
  dev.off()
  log.text(paste(log.description, plot.file, sep=""), log.file=log.file, logging=logging)
}
