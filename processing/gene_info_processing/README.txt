gene_info files are created by parsing Ensembl gtf files for the species.

If a matching genome was available in the list of clusterflow genomes, this was used to get the GC content.
If not, GC content was left blank.

The script used for parsing is create_gene_info_file_from_gtf.pl


output file should contain 
gene_id  gene_name  chromosome  start  end  strand  biotype  length  GC_content  no_of_transcripts

gtf files can be downloaded from ftp://ftp.ensembl.org/pub/current_gtf/....
The script assumes the files are in this format. 

The biotype_families.txt file allows the biotypes in the gtf file to be grouped into families for one of the plots.

to get the biotype families - this isn't species specific - this exact command will work
mysql -uanonymous -P3306 -hensembldb.ensembl.org -Densembl_production_73 -e "select distinct(name),biotype_group from biotype where db_type like '%core%' and is_current=1 order by biotype_group,name;"  > biotype_families.txt