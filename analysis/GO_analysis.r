#rm(list=ls())
#options(encoding="utf-8")
library(data.table)
library(plyr)
# these functions will be packaged up so that the package can just be loaded,
# but for now we'll just source the files
source("/data/private/GOliath/analysis/GOliath_functions.r")
source("/data/private/GOliath/analysis/overrepresentation_test.r")
source("/data/private/GOliath/analysis/plots.r")

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

query.genes <- clean_text(query.genes)
bg.genes <- clean_text(bg.genes)

# file that contains the functional categories and genes within them
species <- as.character(species)
gmt.file.name <- paste0(species, "/", (list.files(species))[1])
print("GO file")
print(gmt.file.name)
gmt.file <- scan(gmt.file.name, sep = "\n", what = "", quiet = TRUE)
print(paste0(length(gmt.file), " categories imported"))

# parse the gmt file
go.categories <- process_GMT(gmt.file)

#===========================
# import the gene info file
#===========================
# this needs sorting properly 
if (grepl(pattern = "Homo_Sapiens", species)) {
    gene_info_file_location <- "/data/private/GOliath/gene_info_data/Homo sapiens/GRCh38.80_gene_info.txt"
    all_gene_info <- fread(gene_info_file_location, select=c(1:5,7:11), stringsAsFactors=TRUE, data.table=FALSE)
} else if (grepl(pattern = "Mus musculus", species)) {
    gene_info_file_location <- "/data/private/GOliath/gene_info_data/Mus musculus/GRCm38.80_gene_info.txt"
    all_gene_info <- fread(gene_info_file_location, select=c(1:5,7:11), stringsAsFactors=TRUE, data.table=FALSE)
} else {
    print("Couldn't find gene info file")
}


#=======================
# for empty background
#=======================
# set the background as all the genes present in the gmt file
# this needs to be a more robust check...
ifelse((length(bg.genes) == 0), bg.genes <- unique(unlist(go.categories)), bg.genes <- remove_duplicates(bg.genes))

#====================
# remove duplicates
#====================
query.genes <- remove_duplicates(query.genes)

#=======================
# convert to upper case
#=======================
#query.genes <- toupper(query.genes)
#bg.genes <- toupper(bg.genes)

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
#go.results <- overrepresentationAnalysis(go.categories, query.genes, bg.genes)

# reduce the number of digits in the output
#go.results$adj.ease.pvalues <- signif(go.results$adj.ease.pvalues, digits=4)

go_results <- overrep_test(go.categories, query.genes, bg.genes)

#==================================
# check against suspect categories
#==================================
#result_table <- read.delim("/data/private/GOliath/jobs/fKnqGgmhXTKGjzSLE5vN/GO_analysis_results.txt")
suspects <- read.delim("/data/private/GOliath/suspect_GO_categories/suspect_categories.txt")

sig_categories <- rownames(go_results)

# in case they categories are in the format "RESPONSE TO CHEMICAL%GOBP%GO:0042221,
# we split by %. If there are no % characters present, it should still work fine.
split_categories <- strsplit(sig_categories, split="%", fixed=TRUE)
ids <- sapply(split_categories, tail, n=1)

# clean up any whitespace so we can do an exact match
flag_locations <- sapply(ids, grep, suspects$GO_ID)

get_description <- function(locations, descriptions){
	paste0(descriptions[locations], collapse=", ")
}

flag_descriptions <- sapply(flag_locations, get_description, suspects$bias_source)

flag_descriptions[flag_descriptions == ""] <- "none found"

go_results$potential_bias <- flag_descriptions

write.table(go_results, "GO_analysis_results.txt", quote=FALSE, sep="\t")

if(!is.null(all_gene_info)){
    
    chromosomes <- list(
        mouse = c(1:19,"MT","X","Y"),
        human = c(1:22,"MT","X","Y"),
        rat = c(1:20,"MT","X","Y"),
        worm = c("I","II","III","IV","MtDNA","V","X"),
        zebrafish = c(1:25,"MT")
    )
    
    # clean up the gene names
    all.gene.info[,"gene_name"] <- toupper(cleanText(all.gene.info[,"gene_name"]))
    
    #=============================================
    # remove any genes not in the background list
    #=============================================
    
    all.gene.info <- all.gene.info[all.gene.info[,"gene_name"] %in% bg.genes,]
    
    # we probably have duplicates in the gene info file but this shouldn't be a problem when we use ensembl ids
    
    # add a TRUE/FALSE column for whether the gene is a query gene
    all.gene.info$query <- all.gene.info[,"gene_name"] %in% query.genes
    
    
    
}


#============
# The plots
#============

#=============
# The GC plot 
#=============

plotting.data <- all.gene.info[,c("GC_content","query")]

png("GC.png")
bar_plot(plotting.data, main = "GC content", xlab = "GC content")
dev.off()














write("", file="finished.flag")

