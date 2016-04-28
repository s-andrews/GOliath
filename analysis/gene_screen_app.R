library(shiny)

#========================
# TODO:
# somehow use genes from a go category as the query genes but not quite sure how to do it....
#   I don't want to mess up the go analysis i.e. it woudl be nice to keep the analysis there
# 
# could we have the option to store multiple sets of query genes - there's no reason why we couldn't add others to the graphs
# though this wouldn't work where we plot the differences in the bar charts
# 
# allow clustering in the main go analysis
# 
# 
# I've had a look at the p values and it all seems to be correct - they just seem quite low compared to other programs (incl giraph)
# I'm not sure whether this has something to do with the background genes used.
# 
# 

#options(shiny.error=browser)

library(data.table)
library(plyr)
library(DT) # for javascript tables
library(VennDiagram)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
#options(shiny.maxRequestSize = 9*1024^2)
options(shiny.maxRequestSize = 18*1024^2)
options(scipen=10000000)


#==========================
# SET UP GLOBAL VARIABLES
#==========================

# any characters that we have for chromosome names - add to as needed
chromosomes <- list(
  mouse = c(1:19,"MT","X","Y"),
  human = c(1:22,"MT","X","Y"),
  rat = c(1:20,"MT","X","Y"),
  worm = c("I","II","III","IV","MtDNA","V","X"),
  zebrafish = c(1:25,"MT")
)

# dataframe of possible species - add to as needed
species.df <- data.frame(common_name= c("mouse", "human", "rat", "worm", "zebrafish"), 
                         scientific_name = c("Mus_musculus", "homo_sapiens", "Rattus", "Caenorhabditis_elegans", "Danio_rerio"))

# species available - this is blank until we check for files that match the species later on 
species.choices <- list("No species info found"="")
  
background.file.choices <- list("")

 # I'm not sure that we need result
go.reactive.values <- reactiveValues(result = 0, filt.categories = 0)


#============
# FUNCTIONS
#============

#===================================================
# check whether we can find a match for the species
#===================================================
speciesInfoFound <- function(textString) if(length(textString) > 1)  sapply(textString, grepText) else grepText(textString)

#===================================================================
# function used within speciesInfoFound function to check whether 
# a match can be found in list.files() for the inputted textString
#===================================================================
grepText <- function(textString) if(sum(grepl(textString, list.files(), ignore.case=T)) > 0) TRUE else FALSE


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
overrepresentationAnalysis <- function(list.of.categories, query.genes, bg.genes, min.query.genes.in.category=3, adj.p.value.threshold=0.05){
  
  query.genes <- toupper(query.genes)
  bg.genes <- toupper(bg.genes)
  
  # check that all the query genes are in the background genes
  query.genes.filt <- query.genes[query.genes %in% bg.genes]
  # filter the list of categories by matching query genes
  # set a min number of query genes for each category - no point having a category with 1 or 2 genes in it.
  categories <- list.of.categories[sapply(list.of.categories, function(x){sum(query.genes.filt %in% x) >= min.query.genes.in.category})]
  
  #=============================================================
  # construct data frame that contains contingency table values
  #=============================================================
  # number of query genes in each category
  df <- data.frame(query.count=sapply(categories, function(x){sum(query.genes.filt %in% x)}))
  # number of background genes in each category
  df$bg.count <- sapply(categories, function(x){sum(x %in% bg.genes)})
  # total number of genes in category
  df$category.lengths <- sapply(categories, length)
  # number of query genes not in category
  query.not.in <- length(query.genes.filt) - df$query.count
  # number of background genes not in category
  bg.not.in <- length(bg.genes) - df$bg.count
  
  # if we wanted to use plain fisher test instead of the ease values
  #df$fisher.pval <- apply(df, 1, function(x){fisher.test(matrix(x,nrow=2))$p.value})
  #df$adj.fisher.pval <- p.adjust(df$fisher.pval,method="BH")
  
  #========================
  # do fisher's exact test 
  #========================
  # use the DAVID stats where they subtract 1 from the number of query genes in the category
  ease.contingency.values <- cbind(df$query.count-1,df$bg.count,query.not.in,bg.not.in)
  
  ease.pvalues <- apply(ease.contingency.values, 1, function(x)fisher.test(matrix(x,nrow=2), alternative = "greater")$p.value)
  df$adj.ease.pvalues <- p.adjust(ease.pvalues, method="BH")
  
  df <- df[df$adj.ease.pvalues <= adj.p.value.threshold,]
  
  df.ordered <- df[order(df$adj.ease.pvalues),]
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

#=======================
# javascript for slider
#=======================
# copied and amended this javascript from http://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
# I don't fully understand it but it seems to work
# There doesn't appear to be a simple way of getting a slider to be on a log scale 
JScode <-
  "$(function() {
setTimeout(function(){
var vals = [1];
var powStart = 4;
var powStop = 8;
for (i = powStart; i <= powStop; i++) {
var val = Math.pow(10, i);
vals.push(val/2000);
vals.push(val/1000);
}
$('#gene_distance').data('ionRangeSlider').update({'values':vals})
}, 5)})"



#==================
# SET UP THE PAGE
#==================

barplot.radio.button.options <- list("plot side-by-side"="raw", "plot differences"="correct")


ui <- fluidPage(

  titlePanel("Gene screen"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(inputId="species", label = "Select species", choices = species.choices),
      
      selectInput(inputId = "background_file_selection", label="Select file containing background gene information", choices=background.file.choices),
      
      textOutput(outputId = "background_file"),
      br(),
	    fileInput(inputId='background_gene_file', label='Upload custom background genes',
	              accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
      br(),
      textInput(inputId = "query_genes", label = "Enter query genes"),     
      actionButton(inputId = "go", label="Go"),
      br(),
	    br(),
      textOutput(outputId = "no_of_genes"),
      br(),
      textOutput(outputId = "duplicate_genes"),
      br(),
      textOutput(outputId = "unmatched_genes")
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel("Gene screen", 
      
          plotOutput(outputId = "GC_density_plot", height=500),
          
          plotOutput(outputId = "length_plot", height=500),
    	  
    	    plotOutput(outputId = "no_of_transcripts_plot", height=500),
          radioButtons(inputId = "transcripts_view", choices=barplot.radio.button.options, inline = TRUE, label=NULL),
    	  
          plotOutput(outputId = "biotype_family_plot", height=500),
    	    radioButtons(inputId = "biotype_family_view", choices=barplot.radio.button.options, inline = TRUE, label=NULL),
          
          plotOutput(outputId = "biotype_plot", height=700),
    	    radioButtons(inputId = "biotype_view", choices=barplot.radio.button.options, inline = TRUE, label=NULL),
          
    	    plotOutput(outputId = "chr_plot", height=500),
    	    radioButtons(inputId = "chr_view", choices=barplot.radio.button.options, inline = TRUE, label=NULL),
          
          plotOutput(outputId = "location_plot", height=500),
          br(),
          tags$head(tags$script(HTML(JScode))),
          # add the tooltip
          tags$div(title="Adjust the minimum distance to closest neighbour", 
                   sliderInput("gene_distance", label= "distance", min = 1, max = 1000000, value = 100))
        ), 
        
        tabPanel("Functional enrichment analysis",
          
          fileInput(inputId="functional_category_file",label="Select file containing functional categories"),
          br(),
          fluidRow(
            column(5,  DT::dataTableOutput("data_table_output"))
          ),  
          br(),
          br(),
          DT::dataTableOutput("overlap_table"),
          textOutput(outputId = "query_genes_in_category"),
    
          plotOutput(outputId = "venn", width="50%")
   
          #plotOutput()
        )
      )  
    )
  )
)

#==================
# THE SERVER CODE
#==================
server <- function(input, output, session) {
 
  # get the species that are available
  species.found <- as.vector(species.df$common_name[rowSums(apply(species.df, 2, speciesInfoFound)) > 0])
  
  # add names to vector
  names(species.found) <- species.found
  
  # update the drop down box with available species
  updateSelectInput(session, inputId="species", choices=as.list(species.found))

  #=================================================
  # we have all.gene.info and background gene info so that all the gene info can go in all.gene.info
  # and if custom background genes are entered then they go into background.gene.info
  # 
  # 
  
  #===========================================================================
  # contains info on location, length, GC content etc for all genes in genome
  # this is only used by background gene info
  #===========================================================================
  all.gene.info <- reactive({
    
    if (nchar(input$background_file_selection) > 1){

      # using fread instead of read.delim as it's so much faster
	    all.gene.info <- fread(input$background_file_selection, select=c(2:5,7:11), stringsAsFactors=TRUE, data.table=FALSE)
      
      all.gene.info[,"gene_name"] <- cleanText(all.gene.info[,"gene_name"])
    }
    else all.gene.info <- NULL
    return(all.gene.info)
  })

  #=================================================================
  # contains info on location, length, GC content etc for all genes 
  # filters if custom background genes have been entered
  # has true/false column for presence of query gene 
  #=================================================================
   background.gene.info <- reactive({
      
    if(is.null(all.gene.info())) background.gene.info <- NULL
    
    else{ 
      background.gene.info <- all.gene.info()
      # this still works with blank query and background fields
  	  background.gene.info$query <- toupper(background.gene.info[,"gene_name"]) %in% toupper(query.genes())
  	  background.gene.info$custom.background <- toupper(background.gene.info[,"gene_name"]) %in% toupper(custom.background.genes())
    }
    # filter if custom background genes have been entered
    if(any(background.gene.info$custom.background)) background.gene.info <- background.gene.info[background.gene.info$custom.background,]
    
    return(background.gene.info)
  })
  
  #=========================== 
  # user entered query genes
  #===========================
  query.genes <- reactive({
    
    return(cleanText(strsplit(input$query_genes, split = " ")[[1]], removeEmpty=TRUE))
    
  })
  
  #=============================== 
  # user entered background genes
  #===============================
  custom.background.genes <- reactive({
    
  	inFile <- input$background_gene_file
  	
  	if(is.null(inFile)) return(NULL) else bg.genes <- scan(inFile$datapath, what="character")  
  	
  })

  #================================== 
  # observe species choice selection
  #==================================
  observeEvent(input$species,{

    if(!(input$species == "No species info found" | input$species == "")){

      matching.files <- c(list.files(pattern=input$species, ignore.case=T), list.files(pattern=as.vector(species.df$scientific_name[species.df$common_name == input$species]),ignore.case=T))
      # filter so we're only including gene_info files
      matching.files <- matching.files[grep(pattern="gene_info", matching.files)]
      names(matching.files) <- matching.files
      updateSelectInput(session, inputId="background_file_selection", choices=as.list(matching.files))
    }
  })
  
  #======================================================================================================  
  # when the GO button is pressed, this text is produced, the analysis doesn't rely on the button though
  #======================================================================================================
  observeEvent(input$go, {
    
    output$no_of_genes <- renderText({
      query.count <- sum(background.gene.info()[,"query"])
      paste(query.count, textPluralCorrection("query genes found in background", "genes", "gene", query.count), "genes")#,  paste(rv$query.genes, collapse = ", "))
	  
    })
    
    output$duplicate_genes <- renderText({
      count <- sum(duplicated(query.genes()))
      paste(count, textPluralCorrection("duplicate genes found in query list:", "genes", "gene", count), paste(query.genes()[duplicated(query.genes())],collapse = ", "))
    })
   
    output$unmatched_genes <- renderText({
      unmatched.genes <- query.genes()[!(toupper(query.genes()) %in% toupper(background.gene.info()[,"gene_name"]))]
      paste(length(unmatched.genes), textPluralCorrection("unmatched genes:", "genes", "gene", length(unmatched.genes)), paste(unmatched.genes, collapse = ", "))
      
    })
  }) 
  
  #================================
  # INFO ABOUT BACKGROUND FILE
  #================================
  output$background_file <- renderText({
    
    if(is.null(all.gene.info())) paste("processing file", input$background_file_selection)
    else	paste(nrow(all.gene.info()), "genes imported from file")
  })
  
    
  #####################
  #
  # GO ANALYSIS PANEL
	#
  ####################
  
  #==================================
  # GO ANALYSIS TABLE 
  #==================================
  observeEvent(input$functional_category_file, {
    
    if(nchar(input$functional_category_file$datapath)>0) {
    
      gmt.file <- scan(input$functional_category_file$datapath, sep="\n", what="")
     # browser()
      categories <- processGMTFile(gmt.file)
      
      result <- overrepresentationAnalysis(categories, query.genes(), background.gene.info()[,"gene_name"])
      
      go.reactive.values$result <- result
      
      go.reactive.values$filt.categories <- categories[rownames(result)]
      
      # the functional enrichment results
      output$data_table_output <- DT::renderDataTable({
        
        result
        
      })      
    }
  })

  #======================================================== 
  # GET THE SELECTED CATEGORIES FROM THE GO ANALYSIS TABLE
  # this is just a reactive value
  # I think the issue here is that the reactive values weren't being accessed because they were within the if statement
  #========================================================
  subset.query <- reactive({

    selected.rows <- input$data_table_output_rows_selected
    
    if(length(go.reactive.values$filt.categories) <= 1 | is.null(input$data_table_output_rows_selected)) return (NULL) else{
      
      subset <-  go.reactive.values$filt.categories[input$data_table_output_rows_selected]
  
      if(length(subset) < 2) return(NULL) else return(sapply(subset, function(x){toupper(query.genes())[toupper(query.genes()) %in% x]}))
    }
 })
 
  #=================================================================
  # TABLE OF SELECTED ROWS FROM GO ANALYSIS TABLE
  # CONTAINS INFO ABOUT NO OF OVERLAPPING GENES BETWEEN CATEGORIES
  #=================================================================
 output$overlap_table <- DT::renderDataTable({
   
   
   if(length(go.reactive.values$filt.categories) > 1){

    if(!(is.null(subset.query()))){

      x <- rbind("total query genes in category" = sapply(subset.query(), length), sapply(subset.query(), function(x){sapply(subset.query(), function(y){sum(x%in%y)})}))
    }
   }
 })
  
  #=========================================================================================================
  # Venn diagram for when 2 categories are selected. If more than 2, overlaps can be looked at in the table	
  #=========================================================================================================
  output$venn <- renderPlot({
   
   print(subset.query())
   
   if(!is.null(subset.query()) & length(subset.query()) == 2){
      # get the number of query genes in each category and the overlap
      nos <- c(length(subset.query()[[1]]), length(subset.query()[[2]]), sum(subset.query()[[1]]%in%subset.query()[[2]]))
      # get the category names
      cat.names <- sapply(strsplit(names(subset.query()),"%"), `[[`, 1)
      # construct the plot
      venn.plot <- draw.pairwise.venn(nos[1], nos[2], nos[3], cat.names, cat.pos=c(0,0), fill=c("#7570B3","#35978F"),lwd=c(0.5,0.5))
      
      grid.draw(venn.plot)
   }
 })
 
  
  #############
  #
  # MAIN PANEL
	#
  #############
  
  #=====================
  # The GC content plot	
  #=====================
  output$GC_density_plot <- renderPlot({
	
	if(is.null(background.gene.info())){
		plot(0:10, type="n", axes = F, ann=F)
		text(x=5, y=5, labels="processing..........", cex=2)
	}
	
    if (length(background.gene.info()) > 1){
      
      plotting.data <- background.gene.info()[,c("GC_content","query")]

      myDensityPlot(plotting.data, "GC content", xlab="GC content")
    }
  })

  #======================
  # The gene length plot	
  #======================
  output$length_plot <- renderPlot({

    if (length(background.gene.info()) > 1){

      plotting.data <- background.gene.info()[,c("length","query")]
      
      myDensityPlot(plotting.data, main="gene length", xlab="log2 length", log=TRUE)
    }  
  })
 
  #============================
  # The no of transcripts plot	
  #============================
  output$no_of_transcripts_plot <- renderPlot({
  
	if (length(background.gene.info()) > 1){   
 
      plotting.data <- background.gene.info()[,c("no_of_transcripts","query")]
      
	   if(input$transcripts_view == "raw") plot.differences = FALSE else plot.differences=TRUE

	   myPercentageBarplot(plotting.data, "no of transcripts per gene", order.numerically = TRUE, plotDifferences = plot.differences)
	}
  })	
 
  #==========================
  # The biotype family plot	
  #==========================
  output$biotype_family_plot <- renderPlot({
    
    if (length(background.gene.info()) > 1){

      plotting.data <- background.gene.info()[,c("biotype_family","query")]
	  
      if(input$biotype_family_view == "raw") plot.differences = FALSE else plot.differences=TRUE
      
	  myPercentageBarplot(plotting.data, "biotype family", plotDifferences = plot.differences)
    }  
  }) 
   
  #===================
  # The biotype plot	
  #=================== 
  output$biotype_plot <- renderPlot({
    
    if (length(background.gene.info()) > 1){
	  
      plotting.data <- background.gene.info()[,c("biotype","query")]
      
      if(input$biotype_view == "raw") plot.differences = FALSE else plot.differences=TRUE
      
	  # some of the biotype categories have very long names but it'd be nice to decide this dynamically
      par(mar =c(20,4,4,2)+0.1)
      myPercentageBarplot(plotting.data, "biotype", las=2, plotDifferences = plot.differences)
    }  
  })
  
  #=====================
  # The chromosome plot	
  #=====================
  # barplot showing number of background and query genes on each chromosome
  output$chr_plot <- renderPlot({
    
    if (length(background.gene.info()) > 1){
     
      plotting.data <- background.gene.info()[,c("chromosome","query")]
      
      if(input$chr_view == "raw") plot.differences = FALSE else plot.differences=TRUE
      
      myPercentageBarplot(plotting.data, "proportion of genes on each chromosome", ordered.categories= chromosomes[[input$species]], plotDifferences = plot.differences)
       
    }    
  })
  
  #=====================
  # Gene proximity plot	
  #=====================
  # barplot showing number of query genes on each chromosome closer than specified distance
  # it only compares the midpoints of genes at the moment
   output$location_plot <- renderPlot({
    
    if(length(query.genes()) < 1){
      plot(0:10, type="n", axes = F, ann=F)
      text(x=5, y=5, labels="enter query genes to display this plot..........", cex=1.5, col="grey")
    }
    
    else{
      
      query.location.data <- background.gene.info()[background.gene.info()$query, c("chromosome", "start", "end")]
      min.distances <- getMinimumDistances(query.location.data)
      
      # convert the slider value
      if(input$gene_distance%%2==0) selected.distance <- (10^(input$gene_distance/2))*1000
      else selected.distance <- ((10^((input$gene_distance+1)/2))/2) *1000
      
      chromosome.counts <- sapply(min.distances, function(x)sum(x<selected.distance))
      
      chromosome.counts.reordered <- chromosome.counts[chromosomes[[input$species]]]

      main.label <- paste("closest neighbour < ", (selected.distance/1000), "kb")
      
      barplot(chromosome.counts.reordered, xlab="chromosome", ylab="no of genes", main=main.label)  
      
    }  
  })
}

shinyApp(ui = ui, server = server)#, session=sessionInfo())





























