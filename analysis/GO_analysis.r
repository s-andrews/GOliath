suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
library(runGOA)

# these functions will be packaged up so that the package can just be loaded,
# but for now we'll just source the files
#source(here::here("analysis","GOliath_functions.r"))
#source(here::here("analysis","overrepresentation_test.r"))
source(here::here("analysis","plots.R"))
source(here::here("analysis","utilities.R"))

# pass in the job folder as the argument
cmdArgs <- commandArgs()

folder.path <- cmdArgs[6]

setwd(folder.path)

# import the config file
config.info <- read.delim("config.txt", header = FALSE, row.names = 1)

type <- config.info["type",]
species <- config.info["species",]

if (is.na(type)) {
    print("gene list type not detected")
}	else if (type == "Unordered") {
    print("Using unordered gene list")	
} 	else if (type == "Ranked") {
    print("Using ranked gene list")
}	else print("Gene list type not recognised")

species_tail <- as.vector(str_split(species, "/", simplify = TRUE)) %>% tail(n=2)
print(paste("Using species", species_tail))

# import the query genes
query_genes <- scan("gene_list.txt", what = "character", quiet = TRUE)
print(paste0(length(query_genes), " query genes imported"))

# import the background genes
bg_genes <- scan("background_list.txt", what = "character", quiet = TRUE)
print(paste0(length(bg_genes), " background genes imported"))

# clean_text removes spaces, characters and converts to upper case
query_genes <- runGOA::clean_text(query_genes)
bg_genes <- runGOA::clean_text(bg_genes)
print(head(query_genes))
print(head(bg_genes))

# file that contains the functional categories and genes within them
#species <- as.character(species)
#gmt.file.name <- paste0(species, "/", (list.files(species))[1])
#print("GO file")
#print(gmt.file.name)
#gmt.file <- scan(gmt.file.name, sep = "\n", what = "", quiet = TRUE)
#print(paste0(length(gmt.file), " categories imported"))
#
## parse the gmt file
#go.categories <- process_GMT(gmt.file)


if (grepl(pattern = "Homo_Sapiens", species)) {
    go.categories <- human_categories
} else if (grepl(pattern = "Mus musculus", species)) {
    go.categories <- mouse_categories
} else {
    print("Couldn't find GO category file")
}



#===========================
# import the gene info file
#===========================
# this needs sorting properly 
if (grepl(pattern = "Homo_Sapiens", species)) {
    gene_info_file_location <- here::here("gene_info_data/Homo sapiens/","GRCh38.80_gene_info.txt")
    all_gene_info <- fread(gene_info_file_location, select = c(1:5,7:11), stringsAsFactors = TRUE, data.table = FALSE)
} else if (grepl(pattern = "Mus musculus", species)) {
    gene_info_file_location <- here::here("gene_info_data/Mus musculus","GRCm38.80_gene_info.txt")
    all_gene_info <- fread(gene_info_file_location, select = c(1:5,7:11), stringsAsFactors = TRUE, data.table = FALSE)
} else {
    print("Couldn't find gene info file")
}

# make sure gene names are upper case
all_gene_info[,"gene_name"] <- toupper(all_gene_info[,"gene_name"])


#=======================
# for empty background
#=======================
# set the background as all the genes present in the gmt file
# this needs to be a more robust check...
ifelse((length(bg_genes) == 0), bg_genes <- unique(unlist(go.categories)), bg_genes <- remove_duplicates(bg_genes))

#====================
# remove duplicates
#====================
query_genes <- remove_duplicates(query_genes)


# check whether all the query genes are in the background genes
query_filt <- query_genes
if (sum(!query_genes %in% bg_genes > 0)) {
    print("not all query genes found in background genes")
    print(query_genes[!query_genes %in% bg_genes])
    
    query_filt <- query_genes[query_genes %in% bg_genes]
}	

#=============================================
# filter options that aren't implemented here
#=============================================
#min.genes.in.category <- 3
#max.genes.in.category <- 5000

go_results <- overrep_test(go.categories, query_filt, bg_genes)

if(is.null(go_results)){
	warning("no significant results found")
	write.table(go_results, "GO_analysis_results.txt", quote = FALSE, sep = "\t")
} else {	

	#==================================
	# check against suspect categories
	#==================================
	suspects <- read.delim(here::here("suspect_GO_categories","suspect_categories.txt"))

	print("read the suspect file")

	sig_categories <- rownames(go_results)

	# it doesn't work just taking the last category as the wikipathways have the organism as the 4th field.
	# get_id <- function(description){
		# ifelse(grepl(x= description, pattern = "%GO:"),
			# strsplit(description, split = "%", fixed = TRUE)[[1]][3], 
			# description
		# )
	# }

	get_id <- function(description){
		if(grepl(x = description, pattern = "%GO:")) {
			strsplit(description, split = "%", fixed = TRUE)[[1]][3] 
		} else if (grepl(x = description, pattern = "reactome", ignore.case = TRUE)){
			strsplit(description, split = "%", fixed = TRUE)[[1]][3]
		} else{
			description
		}	
	}

	ids <- sapply(sig_categories, get_id)
			
	# clean up any whitespace so we can do an exact match
	flag_locations <- lapply(ids, grep, suspects$GO_ID)

	get_description <- function(locations, descriptions){
		paste0(descriptions[locations], collapse = ", ")
	}

	flag_descriptions <- sapply(flag_locations, get_description, suspects$bias_source)

	flag_descriptions[flag_descriptions == ""] <- "none found"

	go_results$potential_bias <- flag_descriptions

	write.table(go_results, "GO_analysis_results.txt", quote = FALSE, sep = "\t")

	#=================
	# summary stats
	#=================
	total_sig_categories    <- nrow(go_results)
	categories_not_flagged  <- sum(flag_descriptions == "none found")
	categories_flagged      <- sum(flag_descriptions != "none found")

  print("printing the size of categories flagged and unflagged")
	print(total_sig_categories)  
	print(categories_not_flagged)
	print(categories_flagged)

	size_of_bias_categories <- suspects %>%
		count(bias_source) %>%
		arrange(desc(n))
		
	#number of GO categories flagged with each bias
	size_of_bias_categories$number_flagged <- sapply(size_of_bias_categories$bias_source, function(y) {
		grep_text <- paste0("\\b", y, "\\b")
		length(grep(x = flag_descriptions, pattern = grep_text))
	})

  df_summary <- data.frame("number of significant categories identified" = total_sig_categories,
             "number of categories flagged as potential biases" = categories_flagged,
             "number of categories not flagged" = categories_not_flagged)


	sink("summary_stats.txt")

# The df_summary would probably be better as some sprintf statements

  # print(paste0("number of significant categories identified = ", total_sig_categories))
  # print(paste0("number of categories flagged as potential biases = ", categories_flagged))       
  # print(paste0("number of categories not flagged = ", categories_not_flagged))
  print(as.data.frame(size_of_bias_categories))

	sink()

  print("I'd like this df info")
  df_summary
  print("And this category size info")
  as.data.frame(size_of_bias_categories)
}

#=================
# screening plots
#=================

if (!is.null(all_gene_info)) {
    
    chromosomes <- list(
        mouse = c(1:19,"MT","X","Y"),
        human = c(1:22,"MT","X","Y"),
        rat = c(1:20,"MT","X","Y"),  
        worm = c("I","II","III","IV","MtDNA","V","X"),
        zebrafish = c(1:25,"MT")
    )
    
    # clean up the gene names
    #   all_gene_info[,"gene_name"] <- toupper(clean_text(all_gene_info[,"gene_name"]))
    
    print(head(all_gene_info))
    
    #=============================================
    # remove any genes not in the background list
    #=============================================
    print("number of background genes =")
    print(length(bg_genes))
    gene_info <- all_gene_info[all_gene_info[,"gene_name"] %in% bg_genes,]
    
    print("number of genes in gene info file that matched the background genes = ")
    print(nrow(all_gene_info))
    
    # we probably have duplicates in the gene info file but this shouldn't be a problem if we use ensembl ids
    
    # add a TRUE/FALSE column for whether the gene is a query gene
    gene_info$query <- gene_info[,"gene_name"] %in% query_filt
    
}


#============
# The plots
#============

#=============
# The GC plot 
#=============

print(head(all_gene_info))

query_GC <- get_GC(query_filt, gene_info)
bg_GC    <- get_GC(bg_genes, gene_info)

print(head(query_GC))

my_plotting_data <- list(query = query_GC, background = bg_GC)

png("GC.png")
density_plot(my_plotting_data, main = "GC content", xlab = "GC content")
dev.off()

#=============
# length plot
#=============
query_lengths <- get_lengths(query_filt, gene_info)
bg_lengths    <- get_lengths(bg_genes, gene_info)

my_plotting_data <- list(query = query_lengths, background = bg_lengths)

png("gene_lengths.png")
density_plot(my_plotting_data, log = TRUE, main = "gene lengths")
dev.off()


#=============
# chr plot
#=============
query_chr <- get_chromosomes(query_filt, gene_info)
bg_chr    <- get_chromosomes(bg_genes, gene_info)

chr_list  <- list(query = query_chr, background = bg_chr)
chr_proportions <- get_chr_percentage(chr_list)

png("chr_plot.png")
bar_plot(chr_proportions, main = "chr")
dev.off()

write("", file = "finished.flag")
