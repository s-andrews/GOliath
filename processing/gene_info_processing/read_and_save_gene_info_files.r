# human
gene_info_file_location <- here::here("gene_info_data/Homo sapiens/","GRCh38.80_gene_info.txt")
human_genfo <- data.table::fread(gene_info_file_location, select = c(1:5,7:11), stringsAsFactors = TRUE, data.table = FALSE)
human_genfo <- mutate(human_genfo, gene_name = toupper(gene_name))
save(human_genfo, file = "processing/gene_info_processing/human_genfo.rda")

# mouse
gene_info_file_location <- here::here("gene_info_data/Mus musculus","GRCm38.80_gene_info.txt")
mouse_genfo <- data.table::fread(gene_info_file_location, select = c(1:5,7:11), stringsAsFactors = TRUE, data.table = FALSE)
mouse_genfo <- mutate(mouse_genfo, gene_name = toupper(gene_name))

save(mouse_genfo, file = "processing/gene_info_processing/mouse_genfo.rda")
