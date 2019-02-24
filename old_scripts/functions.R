library(nlme)

compute.anova <- function(data, response, group.col){
  
  formula.one <- formula(paste0(response," ~ ", group.col))
  response.aov <- summary(aov(formula.one, data=data))
  
  result.row <- data.frame("Metabolite"=response, "ANOVA_PVALUE"=response.aov[[1]][[5]][[1]])

  result.row
}  

# Normalization, value range from 0 to 1
normalization <- function(data, responses){
  
  for (r in responses) {
    
    # values converted to between 0 and 1
    data[,r] <- (data[,r]-min(data[,r]))/(max(data[,r])-min(data[,r]))
  }
  data
}

summarize.results <- function(results, caption="", p.col, p.adj.col, log.file, logging) {
  
  #print(results)
  #p.col <- grep("_p.value$", colnames(results))
  #p.adj.col <-grep("_p.value_adjusted_BH$", colnames(results))
  
  log.text(paste("Summarized results for: ", caption, sep=""), log.file=log.file, logging=logging)
  log.text(paste("- Result rows: ", nrow(results), sep=""), log.file=log.file, logging=logging)
  log.text(paste("- Min p-value: ", min(results[,p.col]), sep=""), log.file=log.file, logging=logging)
  log.text(paste("- Min FDR p-value: ", min(results[,p.adj.col]), sep=""), log.file=log.file, logging=logging)
  log.text(paste("- FDR P < 0.05: ", length(which(results[,p.adj.col] < 0.05)),"\n", sep=""), log.file=log.file, logging=logging)
  
}

# Time-series analysis, compute AUC (area under curve) by subject
perform.auc <- function(data, responses, subject.col, group.col, measurementpoint.col, response.col.title){
  
  require(doMC)
  require(PK)
  
  subjects <- unique(data[, subject.col])
  # Subset individual data
  results <- foreach(h=1:length(subjects), .combine=rbind) %dopar%{
    
    individual <- subjects[h]
    subject.data <- data[which(data[,subject.col] %in% individual),]
    subject.data <- subject.data[,c(subject.col, measurementpoint.col, group.col, responses)]
  

    groups <- unlist(unique(data[,group.col]))
    # Subset data further by group
    group.results <- foreach(j=1:length(groups), .combine=rbind) %dopar% {
      
      group <- groups[j]
      subject.group.data <- subject.data[which(subject.data[,group.col] %in% group),]
      
      # Subset data further by response
      response.result <- foreach(i=1:length(responses), .combine=cbind) %do% {
        
        response <- responses[i]
        response.auc <- auc(conc=subject.group.data[,response], time=subject.group.data[,measurementpoint.col], method='z', design='complete')
        
        # Check
        #print(subject.group.data[,c(measurementpoint.col, response)])
        #print(response.auc[1])
        
        result.col <- data.frame(response.auc[1])
        colnames(result.col) <- response
        
        #if(response.auc > 0.5 & response.auc < 1.0){ print(response.auc) }
          
        result.col
      }
      # Stitch subject and group infront of response row
      temp.df <- data.frame(individual, group)
      colnames(temp.df) <- c(subject.col, "Group")
      response.result <- cbind(temp.df, response.result)
      
      response.result
    }
    group.results
  }
  results
}

compute.anova <- function(data, response, group.col){
  
  formula.one <- formula(paste0(response," ~ ", group.col))
  response.aov <- summary(aov(formula.one, data=data))
  
  result.row <- data.frame("Metabolite"=response, "ANOVA_PVALUE"=response.aov[[1]][[5]][[1]])

  result.row
}  

# Normalization, value range from 0 to 1
normalization <- function(data, responses){
  
  for (r in responses) {
    
    # values converted to between 0 and 1
    data[,r] <- (data[,r]-min(data[,r]))/(max(data[,r])-min(data[,r]))
  }
  data
}

summarize.results <- function(results, caption="", p.col, p.adj.col, log.file, logging) {
  
  #print(results)
  #p.col <- grep("_p.value$", colnames(results))
  #p.adj.col <-grep("_p.value_adjusted_BH$", colnames(results))
  
  log.text(paste("Summarized results for: ", caption, sep=""), log.file=log.file, logging=logging)
  log.text(paste("- Result rows: ", nrow(results), sep=""), log.file=log.file, logging=logging)
  log.text(paste("- Min p-value: ", min(results[,p.col]), sep=""), log.file=log.file, logging=logging)
  log.text(paste("- Min FDR p-value: ", min(results[,p.adj.col]), sep=""), log.file=log.file, logging=logging)
  log.text(paste("- FDR P < 0.05: ", length(which(results[,p.adj.col] < 0.05)),"\n", sep=""), log.file=log.file, logging=logging)
  
}

# Time-series analysis, compute AUC (area under curve) by subject
perform.auc <- function(data, responses, subject.col, group.col, measurementpoint.col, response.col.title){
  
  require(doMC)
  require(PK)
  
  subjects <- unique(data[, subject.col])
  # Subset individual data
  results <- foreach(h=1:length(subjects), .combine=rbind) %dopar%{
    
    individual <- subjects[h]
    subject.data <- data[which(data[,subject.col] %in% individual),]
    subject.data <- subject.data[,c(subject.col, measurementpoint.col, group.col, responses)]
  

    groups <- unlist(unique(data[,group.col]))
    # Subset data further by group
    group.results <- foreach(j=1:length(groups), .combine=rbind) %dopar% {
      
      group <- groups[j]
      subject.group.data <- subject.data[which(subject.data[,group.col] %in% group),]
      
      # Subset data further by response
      response.result <- foreach(i=1:length(responses), .combine=cbind) %do% {
        
        response <- responses[i]
        response.auc <- auc(conc=subject.group.data[,response], time=subject.group.data[,measurementpoint.col], method='z', design='complete')
        
        # Check
        #print(subject.group.data[,c(measurementpoint.col, response)])
        #print(response.auc[1])
        
        result.col <- data.frame(response.auc[1])
        colnames(result.col) <- response
        
        #if(response.auc > 0.5 & response.auc < 1.0){ print(response.auc) }
          
        result.col
      }
      # Stitch subject and group infront of response row
      temp.df <- data.frame(individual, group)
      colnames(temp.df) <- c(subject.col, "Group")
      response.result <- cbind(temp.df, response.result)
      
      response.result
    }
    group.results
  }
  results
}

exclude.samples <- function(data, dataset, identification.col="File", file.path, log.file=log.file, logging=logging){
  # Example
  # exclude.data <- read.csv(paste(input.path, "/data/exclude.csv", sep=""))
  
  # Read the file
  exclude.data <- read.csv(file.path)
  # get value from exclude_samples -column from row where dataset names match
  samples <- exclude.data[which(exclude.data$dataset == dataset),]$exclude_samples
  # split white space separted entries to list
  samples <- unlist(strsplit(as.character(samples), " ", fixed=TRUE))
  # log
  log.text(paste("Removing samples: ", samples, " from dataset ", dataset, sep=""), log.file=log.file, logging=logging) 
  # filter out data from samples(list) with the defined identification.column
  filtered.data <- data[ ! data[,identification.col] %in% samples,]
  # return filtered data
  filtered.data  
}


# Remove rows containing all NAs..
## 1) from data OR
## 2) from data by subset (== Rows that have only NAs in subset)
## print removed rows
## print removed ids (if defined)
### verbose - print removed lines
rm.all.na.rows <- function(data, dataset, subset=NULL, data.id.col=NULL, verbose=TRUE, log.file, logging){
  # Get rows to be removed
  if(is.null(subset)){
    # Only rows with all values NA, from data
    #na.rows <- which(apply(data, 1, function(x) all(is.na(x))))
    na.rows <- which(rowSums(is.na(data)) == ncol(data))
  } else {
    # Only rows with all values NA, from subset
    #na.rows <- which(apply(subset, 1, function(x) all(is.na(x))))
    na.rows <- which(rowSums(is.na(subset)) == ncol(subset))
  }
  # Amount of rows (to be) removed
  n.na.rows <- length(na.rows)
  # Remove na.rows from data
  result.data <- data[-na.rows,]
  # Print na.rows from data
  if(verbose){
    print("Printing removed (NA containing) rows:")
    print(data[na.rows,])
  }
  # Log
  ## Print IDs that are removed
  if(!is.null(data.id.col)){
    ids <- data[na.rows,data.id.col]

    log.text(paste("Removed ", n.na.rows, " samples from dataset ", dataset, ", with following ", data.id.col, " column IDs: \n", sep=""), log.file=log.file, logging=logging)
    log.text(paste(ids, sep=""), log.file=log.file, logging=logging)
  }
  # Print rows
  log.text(paste("Removed ", n.na.rows, " rows from dataset ", dataset, ", with following rows: \n", sep=""), log.file=log.file, logging=logging) 
  log.text(paste(na.rows, sep=""), log.file=log.file, logging=logging)
  
  result.data
}

only.save.results <- function(results, file, data) {
  #results <- cbind(data.frame(Trait=rownames(results)), results)
  write.csv(results, file=file, row.names=FALSE, na="")
}

save.results <- function(results, file, data, add.stats=TRUE) {
  results <- cbind(data.frame(Trait=rownames(results)), results)
  if (add.stats) {
    results <- add.statistics(results, data)
  }
  write.csv(results, file=file, row.names=FALSE, na="")
}

extract.results <- function(results, extract.patterns) {
  
  index <- NULL
  for (p in extract.patterns) {
    index <- c(index, grep(p, colnames(results)))
  }
  index <- unique(index)
  results <- results[,index]
  results <- droplevels(results)
  results
}

# NEEDS REVIEW
# Count N, Mean and SD for each combination of different groups (grouping) per response
summarize <- function(data, responses, grouping, add.zeroes=FALSE, add.ones=FALSE){
  
  require(doMC)
  
  # Make dataframe from groups -> levels as instructions (rows) to handling data
  levelList <- list()
  for(group in rev(grouping)){ levelList[[group]] <- levels(data[,group]) }
  
  fullFactors <- expand.grid(levelList)
  instructions <- rev(fullFactors)
  
  #print(instructions)
  
  # Go through instruction rows one by one
  # subset data by columns and level
  # compute N, Mean, SD
  
  # Function for instruction row handling, computation(N, MEAN, SD)
  rowfunc <- function(row, data, response){
    base.name <- NULL
    
    for(i in 1:length(row)){
      col <- names(row[i])
      level <- as.character(row[i])
      # New column base name by group(s) and level(s)
      if(is.null(base.name)){ base.name <- paste(col, level, sep="_") 
      } else { base.name <- paste(base.name, col, level, sep="_")}
      
      # Subset data
      data <- data[data[col] == level, ]
      
      # Check if data exist (for this instruction row)
      # if not, break the for loop, and data to NULL
      if(nrow(data) == 0){
        data <- NULL
        break
      }
    }
    
    # Continue if data contains something
    if(!is.null(data)){0
                       # Compute N, MEAN, SD
                       m <- mean(data[,response], na.rm=FALSE)
                       sd <- sd(data[,response], na.rm=FALSE)
                       n <- length(which(!is.na(data[,response])))
                       
                       # Stich together
                       result.part <- data.frame(N=n, Mean=m, SD=sd)
                       
                       if (add.zeroes) {
                         zero <- length(which(data[,response] == 0))
                         result.part <- cbind(result.part, data.frame(Zero.N=zero))
                       }
                       if (add.ones) {
                         zero <- length(which(data[,response] == 1))
                         result.part <- cbind(result.part, data.frame(One.N=zero))
                       }
                       colnames(result.part) <- paste(base.name, colnames(result.part), sep="_")
                       
                       if(result.part[1,1] > 0){
                         #print(result.part)
                         return(result.part)
                       } else { return() }
    } else { return() }
  }
  
  summary.df <- data.frame(NULL)
  
  # parallelize
  summary.df <- foreach(i=1:length(responses), .combine=rbind) %dopar% { 
    
    response <- responses[i]
    
    # Subsetting data per response, applying above function
    #for(response in responses){
    
    # Create dataframe with Trait column and response as value
    response.col <- data.frame(Trait=response)
    # Subset data by grouping(list) and response
    data.response <- data[,c(grouping, response)]
    #print(head(data.response))
    # Go through instructions df row by row
    row.dflist <- apply(instructions, 1, function(x) rowfunc(x, data.response, response))
    # Get rid of null vector/list elements
    row.dflist <- row.dflist[!sapply(row.dflist, is.null)]
    # Bind columns together per response
    row <- do.call(cbind, row.dflist)
    # Add response to COMPOUND column
    row <- cbind(response.col, row)
    
    # Keep binding previous df rows to same df for all responses
    #summary.df <- rbind(summary.df, row)
    row
  }
  
  summary.df
}


save.results <- function(results, file, data, add.stats=FALSE, responses, grouping, add.zeroes=FALSE, add.ones=FALSE, logging=FALSE, log.file) {
  results <- cbind(data.frame(Trait=rownames(results)), results)
  if (add.stats) {
    stats <-  summarize(data=data, responses=responses, grouping=grouping, add.zeroes=add.zeroes, add.ones=add.ones)
    results <- merge(x=results, y=stats, by.x="Trait", by.y="Trait")

  }
  log.text(paste("Saving results ", nrow(results), " rows, ", ncol(results), " cols to: ", file, sep=""), logging=logging, log.file=log.file)
  write.csv(results, file=file, row.names=FALSE, na="")
}


replace.values <- function(data, responses, original.value, target.value) {
  
  data[responses][data[responses] == original.value] <- target.value
  data
}




inverse.normalize <- function(data, responses) {
  
  for (r in responses) {
    data[,r] <- qnorm((rank(data[,r], na.last="keep")-0.5) / sum(!is.na(data[,r])))
  }
  data
}


extract.results <- function(results, extract.patterns) {
  
  index <- NULL
  for (p in extract.patterns) {
    index <- c(index, grep(p, colnames(results)))
  }
  index <- unique(index)
  results <- results[,index]
  results <- droplevels(results)
  results
}


mixed <- function(data, response, formula=NULL, random="~ 1 | id", na.action=na.exclude, p.adjust="BH", extract.result=NULL, ci=0.95) {
  
  # Create random term
  if (is.null(random) | is.null(formula)) {
    stop("Formula and Random effect required.")
  }
  original.formula <- formula
  #results <- data.frame()
  
  n <- length(response)
  result.df <- foreach(r=1:length(response), .combine=rbind) %dopar% {
    
    #for(r in 1:length(response)) {
    result.table <- NULL
    print(paste(r, " / ", n, sep=""))
    tmp.response <- response[r]
    print(tmp.response)
    
    #random.term <- paste("~ 1 | ", random, sep="")
    random.term <- random
    
    #formula <- create.formula(tmp.response, covariates, interaction)
    formula <- paste(tmp.response, original.formula, sep="")
    
    
    cat(formula, fill=TRUE)
    l <- NULL
    tryCatch({
      l <- lme(as.formula(formula), data=data, random=as.formula(random.term), na.action=na.action)
      #l <- lme(as.formula(formula), data=data, random=as.formula(random.term), na.action=na.action)
      #l <- withRestarts(lme(as.formula(formula), data=data, random=as.formula(random.term), na.action=na.action), skipError=function() return(NULL))
    }, error = function(e) print(e$message))
    
    if(is.null(l)) {
      cat("Error with LME\n")  
    } else {
      #print(l)
      print(summary(l))
      
      ttable <- data.frame(summary(l)$tTable)
      # Confidence intervals
      if (!is.null(ci)) {
        cintervals <- NULL
        tryCatch({
          cintervals <- intervals(l, level=ci)$fixed[,c(1,3)]
        }, error = function(e) print(e$message))
        
        if (!is.null(cintervals)) {
          colnames(cintervals) <- paste((ci * 100), "_ci_", colnames(cintervals), sep="")
          ttable <- cbind(ttable, cintervals)
        } else {
          cat("Error with LME confidence intervals\n")
          cintervals <- data.frame(lower = rep(NA, times=nrow(ttable)), upper = rep(NA, times=nrow(ttable)) )
          colnames(cintervals) <- paste((ci * 100), "_ci_", colnames(cintervals), sep="")
          ttable <- cbind(ttable, cintervals)
          
        }
        
      }
      result.table <- ttable[1,]
      colnames(result.table) <- paste(rownames(ttable)[1], colnames(result.table), sep="_")
      for(i in 2:nrow(ttable)) {
        tmp.result.table <- ttable[i,]
        colnames(tmp.result.table) <- paste(rownames(ttable)[i], colnames(tmp.result.table), sep="_")
        result.table <- cbind(result.table, tmp.result.table)
      }
      rownames(result.table) <- tmp.response
      #results <- rbind(results, result.table)
    }
    result.table
  }
  
  results <- result.df
  
  # If p.adjust is set, add adjusted p-value columns
  if(!is.null(p.adjust)) {
    r <- results
    p.columns <- grep("p.value$", colnames(r))
    p.columns
    for(h in 1:length(p.columns)) {
      tmp.p <- r[,p.columns[h]]
      tmp.p <- p.adjust(tmp.p, method="BH")
      r <- cbind(r, tmp.p)
      colnames(r)[ncol(r)] <- paste(colnames(r)[p.columns[h]], "_adjusted_", "BH", sep="")
    }
    results <- r
  }
  
  # Order by column names
  results <- results[,order(colnames(results))]
  
  if (!is.null(extract.result)) {
    results <- extract.results(results=results, extract.result=extract.result)
  }
  
  results
}

na.test <-  function (x) {
  w <- sapply(x, function(x)any(is.na(x)))
  if (any(w)) {
    print("NAs in column/column number:")
    print(which(w))
    #stop(paste("NAs in columns", paste(which(w), collapse=", ")))
  }
}

# Replace NA in a certain column in certain position with values from another column (from that position)
replace.na <- function(data, target.col, namereplace.col, replace.value, add.seq.num=TRUE){
  
  # index in column
  x <- which(is.na(data[,target.col]))
  # values of same index in other column
  y <- data[,namereplace.col][x]
  # are all values equal to replace.value
  if(all(y == replace.value)){
    if(add.seq.num){
      # Add sequenctial number to rid off nonunique ID
      y <- paste(as.character(y), seq(1,length(y)), sep="")
    } else {
      y <- as.character(y)
    }
    
    # then replace target.column values with replace.values
    data[,target.col][x] <- y
    # return data
    data
  } else{
    stop("One or more values were not the same as replace.value in namereplace.col!")
  }
}

select.responses <- function(data, pattern, exclude.pattern=NULL, n=NULL) {
  
  responses <- colnames(data)[grep(pattern=pattern, colnames(data))]
  
  index <- NULL
  for (p in exclude.pattern) {
    index <- c(index, grep(p, responses))
  }
  if (is.numeric(index)) { 
    responses <- responses[-index]
  }
  
  if (is.numeric(n)) {
    if (length(responses) < n ) {
      n <- length(responses)
    }
    responses <- sample(responses, n)
  }  
  responses
}



init.log <- function(log.file=NULL, logging=TRUE) {
  if (logging) {
    if(is.null(file)) {
      stop("Log file is not defined.")
    }
    cat("Logging started.\n")
    write(paste(date(), "\n", sep=""), log.file)
  } else {
    cat("Logging is turned off.\n")
  }
}

# Log text 
log.text <- function(text, log.file=NULL, logging=TRUE) {
  
  cat(paste(text, "\n", sep=""))
  if (logging) {
    if(is.null(log.file)) {
      stop("Log file is not defined.")
    }
    write(text, log.file, append=TRUE)
  }
}


# Load metabolite raw data
load.metabo.data <- function(file, scan.lines=20) {
  
  # Skip comment lines starting with #
  tmp.data <- scan(file=file, what="character", nlines=scan.lines, sep="\n")
  index <- grep("^#", tmp.data)
  rm(tmp.data)
  if (length(index) > 0) {
    skip.lines <- max(index)
  } else {
    skip.lines <- 0
  }
  
  if (skip.lines >= scan.lines) {
    stop("Found only header lines.")
  }
  
  data <- read.csv(file, header=TRUE, skip=skip.lines)
  data
  
}

# Extract metabolite compound raw data from the full raw data set
extract.compound.data <- function(data, cols=NULL) {
  
  
  index <- which(colnames(data) == "Compound")
  if (!length(index) > 0) {
    stop("Did not found \"Compound\" column.")
  }
  
  if (is.null(cols)) {
  index <- c(index, grep(".raw.$", colnames(data)))
  }
  else {
    index <- c(index, cols)
  }
  
    
  compound.data <- data[,index]
  
  # Clean sample names
  colnames(compound.data) <- gsub(".raw.$", "", colnames(compound.data))
  
  # Clean metabolite names
  compound.data$Compound <- gsub("@", "a", compound.data$Compound)
  
  # Replace any non-alphanumeric characters with _
  compound.data$Compound <- gsub("[^[:alnum:]]+", "_", compound.data$Compound)
  #compound.data$Compound <- gsub(" ", "_", compound.data$Compound)
  #compound.data$Compound <- gsub("[:punct:]", "_", compound.data$Compound)
  compound.data$Compound <- paste("COMPOUND_", compound.data$Compound, sep="")
  
  # Transpose data
  compound.data.t <- t(compound.data[,-1])
  colnames(compound.data.t) <- compound.data$Compound
  
  # Add DATAFILE column
  result.df <- cbind(data.frame(DATAFILE=rownames(compound.data.t)), compound.data.t)
  
  result.df
  
}


# Load pheno data
load.pheno.data <- function(file) {
  
  data <- read.csv(file, header=TRUE)
  
  required.columns <- c("DATAFILE", "qTOF_SAMPLE_ID", "GROUP")
  if (!all(required.columns %in% colnames(data))) {
    required.columns.text <- paste(required.columns, collapse=", ")
    stop(paste("Missing required columns (", required.columns.text, ")", sep=""))
  }
  
  if(any(is.na(colnames(data)))) {
    warning("Some column names are NA, indicating empty columns")
  }
  data
}


combine.pheno.and.compound.data <- function(pheno.data, compound.data){
  
  if (!isTRUE(all.equal(sort(pheno.data$DATAFILE), sort(compound.data$DATAFILE)))) {
    stop("Problem with matching pheno data to compound data. Likely problem with missing sample/pheno information.")
  }
  result <- merge(x=pheno.data, y=compound.data, by.x="DATAFILE", by.y="DATAFILE")
  result
  
}

# Select columns from data with the given patterns
select.data <- function(data, pattern, log.file, logging) {
  
  selection <- paste(pattern, collapse=", ")
  rows <- nrow(data)
  cols <- ncol(data)
  index <- NULL
  for (p in pattern) {
    index <- c(index, grep(p, colnames(data)))
  }
  data <- droplevels(data[index])
  
  log.text(paste("Subset data:\nSelect cols with patterns: ", selection, "\nOriginal data rows: ", rows, ", cols: ", cols, "\nAfter subsetting, rows: ", nrow(data), ", cols: ", ncol(data), "\n\n", sep=""), logging=logging, log.file=log.file)  
  data
}

rmRows.with.na.rowSums <- function(data){
  # Save rows which rowSums are na
  na.rows <- which(is.na(rowSums(data)))
  # Save data.frame WITHOUT the na.rows
  new.data <- data[-na.rows,]
  # return new.data
  new.data
}
