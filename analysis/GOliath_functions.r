#' Clean text
#'
#' Remove erroneous text characters.
#'
#' @param text character vector
#' @param chars character vector of characters to remove from text
#' @param remove_empty remove empty elements from the vector.
#' @param remove_dup remove duplicated elements from the vector.
#' @param upper_case convert to upper case
#' @return A character vector
#' @examples
#' clean_text("\thello_world")
#' clean_text("\thello_world", upper_case = FALSE)

clean_text <- function(text,
                       chars = c("\t", "\r", "\n", "\\", "\ ", "\""),
                       remove_empty = TRUE,
                       remove_dup = TRUE,
                       upper_case = TRUE
  ) {
      char_pattern <- paste(chars, collapse = "|")
      char_pattern <- paste0("[", char_pattern, "]")
      cleaned <- gsub(pattern = char_pattern, replacement = "", x = as.character(text))
      if (remove_dup) cleaned <- remove_duplicates(cleaned)
      if (upper_case) cleaned <- toupper(cleaned)
      ifelse(remove_empty, return(cleaned[sapply(cleaned, nchar) > 0]), return(cleaned))
}

#' remove_duplicates
#'
#' remove duplicates from vector or data frame
#'
#' @keywords internal
#'
#' @param x vector or data frame to remove duplicates from
#' @param column If a data frame is passed in, the column that is used
#' for deduplicating (default = 1)
#'
#' @return A character vector or data frame.
#' @examples
#' hello_string <- c("h", "e", "l", "l", "o")
#' hello_string
#' remove_duplicates(hello_string)
#'
#' df <- data.frame(hello_string, 1:5)
#' df
#' remove_duplicates(df)
remove_duplicates <- function(x, column = 1) {
  if (is.data.frame(x) | is.matrix(x)) {
    return(x[!duplicated(x[,column]), ])
  } else if (is.list(x)) {
    stop("x must be a vector, data frame or matrix")
  } else if (is.vector(x)) {
    return(unique(x))
  } else{
    warning("data must be a vector, data frame or matrix")
  }
}

#' Import and process GMT file.
#'
#' Import GMT file containing gene sets.
#'
#' The GMT format is described \href{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{GMT format}.
#' The first column should contain gene set names, the 2nd column contains a brief description and the 3rd column contains a whitespace separated list of genes.
#'
#' @param file file path for gmt file containing gene sets
#' @param min_genes minimum number of genes required in the functional category for the category to be loaded.
#' @param max_genes maximum number of genes required in the functional category for the category to be loaded.
#' @return A named list containing all the gene sets imported from \code{gmt_file_path}
#' @examples
#' file <- "http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/Mouse_GO_AllPathways_no_GO_iea_December_01_2018_symbol.gmt"
#' x <- processGMTFile(file, 5, 5000)
process_GMT <- function(file, min_genes = 3, max_genes = 5000) {
	
	gmt_file <- file

#  gmt_file <- data.table::fread(
 #   file,
 #   sep = "\n",
#    header = FALSE,
 #   data.table = FALSE
#  )[, 1]

  categories <- strsplit(gmt_file, "\t")

  names(categories) <- sapply(categories, `[[`, 1)

  genes_in_categories <- lapply(categories, `[`, -(1:2))

  genes_in_categories <- sapply(genes_in_categories, toupper)

  category_lengths <- sapply(categories, length)

  genes_in_categories[category_lengths >= min_genes & category_lengths <= max_genes]
}



