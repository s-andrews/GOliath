# Producing gene info files

In order to produce plots comparing a query set of genes to a background set, we need gene information such as length, location, GC content etc. 
                
The 'gene info', sometimes .genfo is a custom format developed for use with the scripts used within GOliath. Each row is a gene and the column headings are as follows:
```
gene_id  gene_name  chromosome  start  end  strand  biotype  length  GC_content  no_of_transcripts
```

To create a gene info file, a gtf file for the appropriate species is required. This can be downloaded from Ensembl ftp://ftp.ensembl.org/pub/current_gtf/....

The script create_gene_info_file_from_gtf.pl can be run to convert the gtf file to gene info format. The script is written to use Babraham specific file structure to get GC content. If a matching genome is available in the list of clusterflow genomes, this is used to get the GC content. 
If this is unavailable, the GC content is left blank.           

### Biotype families

The biotype_families.txt file allows the biotypes in the gtf file to be grouped into families for one of the plots.

to get the biotype families - this isn't species specific - this exact command will work
```
mysql -uanonymous -P3306 -hensembldb.ensembl.org -Densembl_production_73 -e "select distinct(name),biotype_group from biotype where db_type like '%core%' and is_current=1 order by biotype_group,name;"  > biotype_families.txt
```
