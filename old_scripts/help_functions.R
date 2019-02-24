
extract.mode.responses <- function(responses,mode){
  
  splitted <- strsplit(mode,split="_") %>% unlist()
  chroma <- paste("^COMPOUND_",splitted[1],sep="")
  ion <- splitted[2]
  
  mode.responses <- responses[grepl(chroma,responses) & grepl(ion, responses)]
  mode.responses
}

do.glm <- function(response, predictors, covariate.string=NULL, data, interaction=NULL) {
  
  
  results <- foreach(i=1:length(predictors), .combine=rbind) %dopar% {
    
    
    if (class(predictors) == "list") {
      predictor <- predictors[[i]]
    } else {
      predictor <- predictors[i]
    }
    
    
    if (is.null(covariate.string)) {
      formula <- paste(response, "~", paste(predictor, sep=" + "), sep=" ")
    } else {
      formula <- paste(response, "~", paste(predictor, sep=" + "), covariate.string, sep=" ")
      
    }
    extract <- predictor
    
    if (!is.null(interaction)) {
      formula <- paste(formula, " + ", predictor, "*", interaction, " + ", interaction, sep="")
      if(class(data[,interaction]) == "factor"){
        extract <- paste(predictor, ":", interaction,levels(data[,interaction])[2], sep="")
      }
      else{
        extract <- paste(predictor, ":", interaction, sep="") 
      }
    }
    print(formula)
    print(extract)
    
    result.row <- NULL
    
    tryCatch({
      l <- speedglm(formula=as.formula(formula), family=binomial(), data=data)
      #print(summary(l))
      result.row <- data.frame(Trait=as.character(paste(predictor, sep="_")), LCI95=exp(confint(l)[extract,1]), UCI95=exp(confint(l)[extract,2]), OR=exp(coef(l)[extract]), P=as.numeric(as.character(coef(summary(l))[extract,"Pr(>|z|)"])), stringsAsFactors=FALSE)
    }, error = function(e) print(e$message))
    
    if (is.null(result.row)) {
      result.row <- data.frame(Trait=as.character(paste(predictor, sep="_")), LCI95=NA, UCI95=NA, OR=NA, P=NA, stringsAsFactors=FALSE)
      
    }
    

    
    result.row
  }
  results$"P-FDR" <- p.adjust(results$P, method="BH")
  results
}

perform.anova <- function(data, responses, group.col){
  
  results <- NULL
  results <- foreach(i=1:length(responses), .combine = rbind) %dopar% {
    response = responses[i]
    
    formula.one <- formula(paste0(response," ~ ", group.col))
    response.aov <- summary(aov(formula.one, data=data))
    
    result.row <- data.frame("Metabolite"=response, "ANOVA_PVALUE"=response.aov[[1]][[5]][[1]])
    
    result.row
  }
  results$"P-FDR" <- p.adjust(results$ANOVA_PVALUE, method="BH")
  results
}
