
The script process_geneontology_org_file.pl can process gene ontology data downloaded from geneontology.org. The ontology categories from here should be up-to-date. This script replaces scripts 1 and 2 that are described below. Script 3 - get_GO_parents.pl or get_GO_parents_and_relatives.pl still needs to be run after process_geneontology_org_file.pl in order to get all the GO parents etc.

If ontology information is being downloaded from ensembl, the information below is still all relevant. 

=====================================================================
To get a file of ontology information for all genes in a genome.
The final gmt file can be used in giraph or for the SeqMonk intensity_difference go_category option:

There are 3 perl scripts to use that are all in go_category_processing
1. get_GO_child_terms_ensembl_api.pl
2. process_go_child_terms.pl
3. get_GO_parents.pl or get_GO_parents_and_relatives.pl

The only part of this that is species specific is the initial step - downloading the GO category information from Ensembl.

See below for details for each script.



1. Use the script get_GO_child_terms_ensembl_api.pl - this does not iterate up and get all the parent terms, that part is commented out in the script.
It connects to the ensembl api and gets all the lowest level terms for each gene.
module load bioperl 
module load ensemblapi

USAGE: perl get_GO_child_terms_ensembl_api.pl Homo_sapiens > output_file.txt

output e.g.
YKL222C GO:0005515      molecular_function      protein binding

2. This needs to be processed using the process_go_child_terms.pl script. the output of this is in gmt format.
    eg. HABITUATION%GOBP%GO:0046959     habituation     cat-2   magi-1  dop-1

perl process_go_child_terms.pl input_file.txt
	
3. To get all the parent GO categories, use get_GO_parents.pl or get_GO_parents_and_relatives.pl
The parents script only includes the is_a relationships.
The parents_and_relatives script also includes the other relationships, i.e. part_of, positively_regulates, negatively_regulates and regulates.
The output of this is a new gmt file that includes the child and parent (relatives).

It relies on the go-basic.obo file which needs to be updated frequently, it is released daily - download from http://purl.obolibrary.org/obo/go/go-basic.obo

perl get_GO_parents.pl input_file.txt
perl get_GO_parents_and_relatives.pl input_file.txt




to check the available databases
mysql -u anonymous -h ensembldb.ensembl.org -P 3306
SHOW DATABASES LIKE "%core%";
