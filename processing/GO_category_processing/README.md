# Processing Gene ontology Information

In order to perform gene set analyses we need up-to-date files containing gene ontology categories and the genes mapped to these categories for a number of species. Ontology information for genes is updated and released fairly regularly from a number of sources.

## GMT - Gene Matrix Transposed file format 

GMT is a file format originially used by the Broad Institute for their Gene Set Enrichment Analysis https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29               

In this format there is a separate line for each ontology category. The first column contains a name, the second a brief description and the third and any subsequent columns contain the gene(s) belonging to that category
```
HABITUATION%GOBP%GO:0046959                    habituation                    cat-2   magi-1   dop-1                   
RESPIRATORY CHAIN COMPLEX I%GOCC%GO:0045271    respiratory chain complex      gas-1   lpd-5    nuo-1   C18E9.4   C33A12.1          
```
### Downloading GMT files directly

Gene set information for human, mouse and rat can be downloaded from http://download.baderlab.org/EM_Genesets/. These are currently released monthly. 

GMT files containing pathway information can be downloaded directly from pathway commons http://www.pathwaycommons.org/archives/

### Producing GMT files

Within this GO_category_processing folder are scripts to process ontology data and convert it to GMT format. This allows more species to be included and files to be updated when required. 

There are multiple places that ontology information can be obtained from.


## Processing data from geneontology.org

#### Gene association files for a number of species can be downloaded from http://geneontology.org/page/download-annotations 

The .gaf.gz file should be downloaded and the the following scripts need to be used in sequence:             
1. process_geneontology_org_file.pl 
2. get_GO_parents.pl or get_GO_parents_and_relatives.pl as described above

(There are other file formats on the web page for some species - not sure whether the script works for these.)

#### 1. process_geneontology_org_file.pl  

The perl script processes the file and converts it to GMT format. It only contains the most specific categories i.e. the lowest level child categories. To populate the file with all the parent categories, another script needs to be run.

Usage:


#### 2. get_GO_parents.pl OR get_GO_parents_and_relatives.pl

The parents script only includes the is_a relationships.
The parents_and_relatives script also includes the other relationships, i.e. part_of, positively_regulates, negatively_regulates and regulates.
The output of this is a new gmt file that includes the child and parent (relatives).

Usage:  perl get_GO_parents.pl input_file.txt
	perl get_GO_parents_and_relatives.pl input_file.txt

A go-basic.obo file is required - by default the script looks for this in the current working directory.

### go-basic.obo
To process the parent and child mappings, the file go-basic.obo is required. This is released daily from http://purl.obolibrary.org/obo/go/go-basic.obo



## Processing data from Ensembl

Ontology information can be downloaded from Ensembl instead of geneontology.org. 
The following scripts need to be used in sequence
1. get_GO_child_terms_ensembl_api.pl 
2. process_go_child_terms.pl
3. get_GO_parents.pl or get_GO_parents_and_relatives.pl as described above

#### 1. Using the Ensembl API - get_GO_child_terms_ensembl_api.pl 

The script get_GO_child_terms_ensembl_api.pl connects to the ensembl api and gets all the lowest level terms for each gene.
module load bioperl 
module load ensemblapi

USAGE: perl get_GO_child_terms_ensembl_api.pl Homo_sapiens > output_file.txt

output e.g.
YKL222C GO:0005515      molecular_function      protein binding

to check the available databases
mysql -u anonymous -h ensembldb.ensembl.org -P 3306
SHOW DATABASES LIKE "%core%";

In theory this script could iterate up and get all the parent terms, but this takes an excessively long times so that part is commented out in the script

#### 2. Process the downloaded information - process_go_child_terms.pl

The output of the api script needs to be processed.

perl process_go_child_terms.pl input_file.txt

#### 3. Run GO_parents.pl or get_GO_parents_and_relatives.pl

The parents script only includes the is_a relationships.
The parents_and_relatives script also includes the other relationships, i.e. part_of, positively_regulates, negatively_regulates and regulates.
The output of this is a new gmt file that includes the child and parent (relatives).

Usage:  perl get_GO_parents.pl input_file.txt
	perl get_GO_parents_and_relatives.pl input_file.txt

A go-basic.obo file is required - by default the script looks for this in the current working directory.

### go-basic.obo
To process the parent and child mappings, the file go-basic.obo is required. This is released daily from http://purl.obolibrary.org/obo/go/go-basic.obo

