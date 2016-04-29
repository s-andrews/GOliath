library(data.table)
library(plyr)

# hardcoded bits that would need changing
#gene.info.filename <- ("/home/bigginsl/GOliath/jobs/test/Mouse_Mus_musculus.GRCm38.80_gene_info.txt")

#this.species <- "mouse" # for correct chromosomal ordering

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
query.genes <- scan("gene_list.txt", what="character")
print(paste0(length(query.genes), " query genes imported"))

# import the background genes
bg.genes <- scan("background_list.txt", what="character")
print(paste0(length(bg.genes), " background genes imported"))

# file that contains the gene properties
#species <- as.character(species)
#gene.info.file.name <- paste0(species, "/", (list.files(species))[2])
#print("gene info file")
#print(gene.info.file.name)

gene.info.filename <- ("/home/bigginsl/GOliath/jobs/test/Mouse_Mus_musculus.GRCm38.80_gene_info.txt")

all.gene.info <- fread(gene.info.filename, select=c(1:5,7:11), stringsAsFactors=TRUE, data.table=FALSE)

print(paste0(nrow(all.gene.info), " genes imported"))

this.species <- "mouse"


#=============================================================
# density plot used for displaying GC content and gene length
#=============================================================
myDensityPlot <- function(plotting.data, main="", xlab="", log=FALSE){
  
  # need at least 2 values to plot density
  if(is.null(plotting.data) | nrow(plotting.data)<2) return
  
  if(all(is.na(plotting.data[,1])) | nrow(plotting.data) < 1){
    plot(0:10, type="n", axes = F, ann=F)
    text(x=5, y=5, labels=paste("no information available for", main))
    return()
  }
  
  # remove any NA values
  plotting.data <- na.omit(plotting.data)
  
  background.values <- plotting.data[,1]
  if(sum(plotting.data[,2])>1) query.values <- plotting.data[,1][plotting.data[,2]] else query.values <- NULL
  
  if(log == TRUE){
    background.values <- log2(background.values)
    
    # we allow null query values (for when a genelist hasn't been loaded) but we'd get an error if we tried to log transform them
    if(!is.null(query.values)) query.values <- log2(query.values)
  }
  
  if(!is.null(query.values))  density.query <- density(query.values, na.rm=T) else density.query <- NULL
  
  x.range <- range(c(query.values, background.values))
  y.range <- range(c(density.query)$y, density(background.values, na.rm=T)$y)
  
  plot(density(background.values, na.rm=T), main=main, xlim=x.range, ylim=y.range, lwd=2, xlab=xlab)
  lines(density.query,  col="red", lwd=2)
  
  legend("topright", legend = c("background", "query"), fill=c("black", "red"), bty="n")  
}

#===============================================
# the barplot used for biotypes and chromosome 
#===============================================
myPercentageBarplot <- function(plotting.data, main="", xlab="", las=1, order.numerically=FALSE, ordered.categories=NULL, plotDifferences=FALSE){
  
  if(is.null(plotting.data)) return()
  
  if(is.null(nrow(plotting.data)) | nrow(plotting.data) < 1){
    plot(0:10, type="n", axes = F, ann=F)
    text(x=5, y=5, labels=paste("no information available for", main), cex=2)
    return()
  }
  
  plot.title <- colnames(plotting.data)[1]
  
  # count function is from plyr package
  bg.counts <- count(plotting.data, plot.title)
  query.counts <- count(plotting.data, plot.title, "query")
  
  count.matrix <- matrix(data=c(bg.counts$freq, query.counts$freq), ncol=2)
  
  rownames(count.matrix) <- bg.counts[,1]
  
  percentage.matrix <- apply(count.matrix, 2, function(x){if(sum(x) == 0) total<-1  else total <- sum(x); (x/total)*100})
  
  if(order.numerically==FALSE & is.null(ordered.categories)) ordered.percentage.matrix <- percentage.matrix
  
  else{
    
    if(order.numerically==TRUE) reordered <- as.character(sort(type.convert(rownames(percentage.matrix)))) 
    else if(!is.null(ordered.categories)) reordered <- ordered.categories[ordered.categories%in%rownames(percentage.matrix)]
    
    ordered.percentage.matrix <- percentage.matrix[reordered,]
  }
  
  
  # if there's only one column we still need to make it into a matrix
  if(is.null(dim(ordered.percentage.matrix)) & length(ordered.percentage.matrix) > 0) {
    ordered.percentage.matrix <- t(ordered.percentage.matrix)
    rownames(ordered.percentage.matrix) <- bg.counts[,1]
  }  
  
  if(plotDifferences == FALSE){
    barplot(t(ordered.percentage.matrix), beside=T, main=main, xlab=xlab, ylab="%", las=las)#,  names.arg = levels(background.values))
    legend("topright", legend = c("background", "query"), fill=c("black", "grey"), bty="n")  
  }
  else{
    diff.data <- ordered.percentage.matrix[,2] - ordered.percentage.matrix[,1]
    barplot(diff.data, main=main, xlab=xlab, ylab="query - background proportion", las=las)
  }
}

#=================================================
# function to get minimum distances between genes 
#=================================================
getMinimumDistances <- function(query.location.data){
  
  # get the number of genes per chromosome
  chr.counts <- table(query.location.data[,"chromosome"])
  
  # remove the chromosomes where there's only 1 gene
  query.location.data <- query.location.data[query.location.data[,"chromosome"] %in% names(chr.counts)[chr.counts > 1],]
  
  centrepoints <- query.location.data[,"start"] + (query.location.data[,"end"]-query.location.data[,"start"])/2
  
  # get the distance between a gene and its closest neighbour
  min.distances <- tapply(centrepoints, INDEX = query.location.data[,"chromosome"], FUN = function(x){
    
    if(length(x)<2) return (NA) else{
      
      my.dist <- as.matrix(dist(x, upper=T))
      # we don't want the diagonals to be 0
      diag(my.dist) <- NA
      
      apply(my.dist, 1, min, na.rm=T)
    }
  })
  return(min.distances)
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

#=======================
# for empty background
#=======================
# set the background as all the genes present in the gene info file
# this needs to be a more robust check...
ifelse((length(bg.genes) == 0), bg.genes <- all.gene.info$gene_name, bg.genes <- removeDuplicates(bg.genes, "background genes"))

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


#=====================================================================
# any characters that we have for chromosome names - add to as needed
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



#============
# The plots
#============

#pdf("/home/bigginsl/GOliath/jobs/test/screening_plots.pdf")

setwd("/home/bigginsl/GOliath/jobs/test/")

#=============
# The GC plot 
#=============

plotting.data <- all.gene.info[,c("GC_content","query")]

png("GC.png")
myDensityPlot(plotting.data, "GC content", xlab="GC content")
dev.off()

#======================
# The gene length plot	
#======================

plotting.data <- all.gene.info[,c("length","query")]

png("gene_lengths.png")
myDensityPlot(plotting.data, main="gene length", xlab="log2 length", log=TRUE)
dev.off()

#=========================
# transcript number plot
#========================
plotting.data <- all.gene.info[,c("no_of_transcripts","query")]
# this was for toggling the view
#if(input$transcripts_view == "raw") plot.differences = FALSE else plot.differences=TRUE

png("transcripts_per_gene.png")
par(mfrow=c(2,1))
myPercentageBarplot(plotting.data, "no of transcripts per gene", order.numerically = TRUE, plotDifferences = FALSE)
myPercentageBarplot(plotting.data, "no of transcripts per gene", order.numerically = TRUE, plotDifferences = TRUE)
dev.off()

#==========================
# The biotype family plot	
#==========================
plotting.data <- all.gene.info[,c("biotype_family","query")]

png("biotype_families.png")
par(mfrow=c(2,1))
myPercentageBarplot(plotting.data, "biotype family", plotDifferences = FALSE)
myPercentageBarplot(plotting.data, "biotype family", plotDifferences = TRUE)
dev.off()

#===================
# The biotype plot	
#=================== 
plotting.data <- all.gene.info[,c("biotype","query")]
# some of the biotype categories have very long names but it'd be nice to decide this dynamically

png("biotypes.png")
par(mfrow=c(2,1))
par(mar =c(10,4,4,2)+0.1)
myPercentageBarplot(plotting.data, "biotype", las=2, plotDifferences = FALSE)
myPercentageBarplot(plotting.data, "biotype", las=2, plotDifferences = TRUE)
dev.off()

#=====================
# The chromosome plot	
#=====================
# barplot showing number of background and query genes on each chromosome

plotting.data <- all.gene.info[,c("chromosome","query")]

png("chromosomes.png")
par(mar =c(5,4,4,2)+0.1)
par(mfrow=c(2,1))
myPercentageBarplot(plotting.data, "proportion of genes on each chromosome", ordered.categories= chromosomes[[this.species]], plotDifferences = FALSE)
myPercentageBarplot(plotting.data, "proportion of genes on each chromosome", ordered.categories= chromosomes[[this.species]], plotDifferences = TRUE)

dev.off()
