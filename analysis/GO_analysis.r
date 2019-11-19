suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(runGOA))
suppressMessages(library("magrittr"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))

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
tibble::as_tibble(t(config.info)) -> tidy_config

type    <- pull(tidy_config, type)
species <- pull(tidy_config, species)
min_category_size <- pull(tidy_config, min_category_size)
max_category_size <- pull(tidy_config, max_category_size)

if(min_category_size == ""){
  min_category_size = 10
  print("Setting minimum category size to 10 as no information found from input form")
}

if(max_category_size == ""){
  max_category_size = 500
  print("Setting maximum category size to 500 as no information found from input form")
}

min_category_size <- as.numeric(min_category_size)
max_category_size <- as.numeric(max_category_size)
 
#type <- config.info["type",]
#species <- config.info["species",]
#min_category_size <- config.info["min_category_size	",]
#min_category_size <- config.info["min_category_size	",]

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
bg_genes    <- runGOA::clean_text(bg_genes)
print(head(query_genes))
print(head(bg_genes))

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


#===========================
# The GO analysis
#=============================
go_results <- runGOA::overrep_test(go.categories, query_filt, bg_genes, min_genes_in_category = min_category_size, max_genes_in_category = max_category_size)

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

gc_data <- tibble(GC = query_GC, type = "query") %>%
  bind_rows(tibble(GC = bg_GC, type = "background"))


p <- ggplot(gc_data, aes(x = GC, fill = type, color = type)) +
  geom_density(alpha = 0.5, size = 1.5) +
  ggtitle("\nGC content of genes\n") +
  scale_color_manual(values = c("black", "red3")) +
  scale_fill_manual(values = c("black", "red3")) +
  theme(
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16),
    axis.title   = element_text(size = 20),
    axis.text    = element_text(size = 14),
    title        = element_text(size = 22),
    legend.spacing.x = unit(0.2, 'cm')
  )

png("GC.png", width = 600, height = 400)
p
dev.off()


#=============
# length plot
#=============
query_lengths <- get_lengths(query_filt, gene_info)
bg_lengths    <- get_lengths(bg_genes, gene_info)

length_data <- tibble(gene_length = query_lengths, type = "query") %>%
  bind_rows(tibble(gene_length = bg_lengths, type = "background"))

p <- ggplot(length_data, aes(x = log2(gene_length), fill = type, color = type)) +
  geom_density(alpha = 0.5, size = 1.5) +
  ggtitle("\nGene lengths\n") +
  scale_color_manual(values = c("black", "red3")) +
  scale_fill_manual(values = c("black", "red3")) +
  labs(x = "log2 gene length") +
  theme(
    legend.title = element_blank(), 
    legend.text  = element_text(size = 16),
    axis.title   = element_text(size = 20),
    axis.text    = element_text(size = 14),
    title        = element_text(size = 22),
    legend.spacing.x = unit(0.2, 'cm')
  )

png("gene_lengths.png", width = 600, height = 400)
p
dev.off()

#=============
# chr plot
#=============
query_chr <- get_chromosomes(query_filt, gene_info)
bg_chr    <- get_chromosomes(bg_genes, gene_info)

chr_list  <- list(query = query_chr, background = bg_chr)

zz <- file("all.Rout", open="wt")
sink(zz, type="message")

chr_proportions <- get_chr_percentage(chr_list)

## reset message sink and close the file connection
sink(type="message")
close(zz)

chr <- rownames(chr_proportions)

tidy_chr <- tibble::as_tibble(chr_proportions) %>%
  mutate(chr = factor(chr, levels = chr)) %>%
  tidyr::gather(`query`, `background`, key = "type", value = "proportion")

p <- ggplot(tidy_chr, aes(x = chr, y = proportion, fill = type, color = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  ggtitle("\nChromosomal locations\n") +
  coord_flip() +
  labs(y = "% of genes on chromosome") +
  scale_color_manual(values = c("black", "red3")) +
  scale_fill_manual(values = c("black", "red3")) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 16),
    axis.title   = element_text(size = 20),
    axis.text    = element_text(size = 14),
    axis.title.y = element_text(angle = 0,  hjust = 0),
    title        = element_text(size = 22),
    legend.spacing.x = unit(0.2, 'cm'))


png("chr_plot.png", width = 600, height = 700)
p
dev.off()

write("", file = "finished.flag")

#warning("Warning message 6")