#rm(list=ls())
#options(encoding="utf-8")

# pass in the job folder as the argument

cmdArgs <- commandArgs()

folder.path <- cmdArgs[6]

setwd(folder.path)

# import the config file
config.info <- read.delim("config.txt", header=FALSE, row.names=1)

type <- config.info["type",]
species <- config.info["species",]

if(is.na(type)){
	print("gene list type not detected")
}	else if(type == "Unordered"){
		print ("Unordered gene list found")	
} 	else if(type == "Ranked"){
		print ("Ranked gene list found")
}	else print("Gene list type not recognised")

print(paste("Using species", species))

# import the query genes
query.genes <- scan("gene_list.txt", what="character", quiet=TRUE)
print(paste0(length(query.genes), " query genes imported"))

# import the background genes
bg.genes <- scan("background_list.txt", what="character", quiet=TRUE)
print(paste0(length(bg.genes), " background genes imported"))

# file that contains the functional categories and genes within them
species <- as.character(species)
gmt.file.name <- paste0(species, "/", (list.files(species))[1])
print("GO file")
print(gmt.file.name)
gmt.file <- scan(gmt.file.name, sep="\n", what="", quiet=TRUE)
print(paste0(length(gmt.file), " categories imported"))


#############
# functions
#############

#===================
# remove duplicates
#===================
removeDuplicates <- function(string, name.of.string){
	if(sum(duplicated(string)) > 0){
		string <- unique(string)
		print(paste("Duplicates found and removed.", length(string), name.of.string, "remaining."))
	}
	return(string)
}

#===================================
# function to process the gmt file
#===================================
processGMTFile <- function(gmt.file, min.genes.in.category=3, max.genes.in.category=5000){
 
  # extract the category info
  categories <- strsplit(gmt.file, "\t")
  # Extract the first vector element and set it as the list element name
  names(categories) <- sapply(categories, `[[`, 1)
  # Remove the first 2 vector elements from each list element
  categories <- lapply(categories, `[`, -(1:2))
  # convert to uppercase so we can match to query and bg genes
  categories <- sapply(categories, toupper)
  # top-level filter - remove very large and small categories
  categories[sapply(categories, length) >= min.genes.in.category & sapply(categories, length) <= max.genes.in.category]
  
}

#=================
# THE GO ANALYSIS
#=================

# list of categories is all the functional categories from the gmt file
overrepresentationAnalysis <- function(list.of.categories, query.genes, bg.genes, min.query.genes.in.category=3, adj.p.value.threshold=0.05, ease=TRUE){
  

  # filter the list of categories by matching query genes
  # set a min number of query genes for each category - no point having a category with 1 or 2 genes in it.
  categories <- list.of.categories[sapply(list.of.categories, function(x){sum(query.genes %in% x) >= min.query.genes.in.category})]
  
  #=============================================================
  # construct data frame that contains contingency table values
  #=============================================================
  # number of query genes in each category
  df <- data.frame(query.count=sapply(categories, function(x){sum(query.genes %in% x)}))
  # number of background genes in each category
  df$bg.count <- sapply(categories, function(x){sum(x %in% bg.genes)})
  # total number of genes in category
  df$category.lengths <- sapply(categories, length)
  # number of query genes not in category
  query.not.in <- length(query.genes) - df$query.count
  # number of background genes not in category
  bg.not.in <- length(bg.genes) - df$bg.count 
  # number of categories that were discarded as contained too few query genes
  untested_categories <- length(list.of.categories) - length(categories)	
	
	
  if(ease == FALSE){	
  #========================
  # do fisher's exact test 
  #========================
    contingency.values <- cbind(df$query.count,df$bg.count,query.not.in,bg.not.in)

    df$fisher.pval <- apply(contingency.values, 1, function(x){fisher.test(matrix(x,nrow=2), alternative = "greater")$p.value})

	   # slight fudge because we need to do the multiple testing correction for more categories than we have here.
	   # we have only actually tested categories with a certain number of query genes in them so as not to slow 
	   # the code down with extra fishers tests that wouldn't be biologically significant

    extended.pvals <- c(df$fisher.pval, rep(10, untested_categories))
    adjusted_extended <- p.adjust(extended.pvals,method="BH")

    # take the top 
    df$adj.fisher.pval <- adjusted_extended[1:length(categories)]
	  
    df <- df[df$adj.fisher.pval <= adj.p.value.threshold,]
    return(df[order(df$adj.fisher.pval),])
  }	  
  
  else{	
  #========================
  # use the DAVID stats where they subtract 1 from the number of query genes in the category
    ease.contingency.values <- cbind(df$query.count-1,df$bg.count,query.not.in,bg.not.in)

    df$ease.pval <- apply(ease.contingency.values, 1, function(x)fisher.test(matrix(x,nrow=2), alternative = "greater")$p.value)

    extended.ease.pvals <- c(df$ease.pval, rep(10, untested_categories))
    adjusted_extended.ease <- p.adjust(extended.ease.pvals,method="BH")

    df$adj.ease.pvalues <- adjusted_extended.ease[1:length(categories)]

    df <- df[df$adj.ease.pvalues <= adj.p.value.threshold,]

    return(df[order(df$adj.ease.pvalues),])
  }				 
}

#=====================================================
# function to strip line endings and other characters
# removeEmpty=T will remove blank items in the vector
#======================================================
cleanText <- function(character.vector, removeEmpty=F){
  
  chars <- "[\t|\r|\n| ]"
  
  cleaned <- gsub(pattern = chars, replacement = "", x=as.character(character.vector))
  
  if(removeEmpty) cleaned[sapply(cleaned, nchar) > 0] else cleaned
}

#==========================================================================================================================================================
# function to correctly pluralise text depending on number eg 1 gene, 2 genes
# write text output within server assuming the number of genes etc is greater than 1, then use this function to check and correct the text if number <= 1 
# e.g. textPluralCorrection(textString="duplicate genes removed:", pluralString="genes", singularString="gene", count=count)
#==========================================================================================================================================================
textPluralCorrection <- function(textString, pluralString, singularString, count){
  
  if(count == 1) gsub(pluralString, singularString, textString) else if(count==0 | is.null(count)) gsub(":$", "", textString) else textString
}

# parse the gmt file
go.categories <- processGMTFile(gmt.file)

#=======================
# for empty background
#=======================
# set the background as all the genes present in the gmt file
# this needs to be a more robust check...
ifelse((length(bg.genes) == 0), bg.genes <- unique(unlist(go.categories)), bg.genes <- removeDuplicates(bg.genes, "background genes"))

#====================
# remove duplicates
#====================
query.genes <- removeDuplicates(query.genes, "query genes")

#=======================
# convert to upper case
#=======================
query.genes <- toupper(query.genes)
bg.genes <- toupper(bg.genes)

# check whether all the query genes are in the background genes
if(sum(!query.genes %in% bg.genes > 0)){
	print("not all query genes found in background genes")
	print(query.genes[!query.genes %in% bg.genes])
	
	query.genes <- query.genes[query.genes %in% bg.genes]
}	
#=============================================
# filter options that aren't implemented here
#=============================================
#min.genes.in.category <- 3
#max.genes.in.category <- 5000

# perform the Fishers Exact test to get results
go.results <- overrepresentationAnalysis(go.categories, query.genes, bg.genes)

# reduce the number of digits in the output
go.results$adj.ease.pvalues <- signif(go.results$adj.ease.pvalues, digits=4)

write.table(go.results, "GO_analysis_results.txt", quote=FALSE, sep="\t")
write("", file="finished.flag")

