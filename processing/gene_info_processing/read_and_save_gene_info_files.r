# human
gene_info_file_location <- here::here("gene_info_data/Homo sapiens/","GRCh38.80_gene_info.txt")
human_genfo <- fread(gene_info_file_location, select = c(1:5,7:11), stringsAsFactors = TRUE, data.table = FALSE)
save(human_genfo, "human_genfo.rda")

# mouse
gene_info_file_location <- here::here("gene_info_data/Mus musculus","GRCm38.80_gene_info.txt")
mouse_genfo <- fread(gene_info_file_location, select = c(1:5,7:11), stringsAsFactors = TRUE, data.table = FALSE)
save(mouse_genfo, "mouse_genfo.rda")
