#!/usr/bin/perl

# processes the gene ontology data downloaded from geneontology.org - http://geneontology.org/page/download-annotations
# input files should be consistent across species but we should check that the information in each column makes sense - yhis is not yet implemented
# This script will use the gene IDs provided in the input file from geneontology.org and mappings to ensembl ids etc can be carried out separately 
# as the process will vary depending on the organism

# It relies on the go-basic.obo file which needs to be updated frequently, 
# it is released daily - download from http://purl.obolibrary.org/obo/go/go-basic.obo

# usage perl process_geneontology_org_file --obo go-basic.obo --gene_ontology_file gene_association.mgi

# input file format eg.
# MGI     MGI:1915571     0610007P14Rik           GO:0005783  

# produces output like the gmt files used at the Broad and Bader lab e.g.
# TRANSPORTER COMPLEX%GOCC%GO:1990351     transporter complex     ccb-2   ccb-1

# To get the GO parents (and relatives) a separate script can be run once this has completed.

use warnings;
use strict;
use Getopt::Long;

my $obo_file;
my $gene_ontology_file;
my %go_categories;
# used to store all the go ids and descriptions from the obo file.
my %go_id_description;
# key is child, values are parents
my %child_parents;

my $config_result = GetOptions(
"obo=s" => \$obo_file,
"gene_ontology_file=s"  => \$gene_ontology_file
#"help" => \$help
);
die "Could not parse options" unless ($config_result);

print "\nobo file is $obo_file \n";
load_obo_file_into_hash($obo_file);


my $outfile = $gene_ontology_file;
$outfile =~ s/.txt//;
$outfile = "$outfile.go_terms.txt";
open (OUT, ">", $outfile) or die  "Couldn't open output file: $! \n";

open (IN, $gene_ontology_file) or die "Cannot open file $gene_ontology_file\n";

while (my $line = <IN>){
	
	next if($line =~ /^!/);
	
	my @fields = split(/\t/, $line);
	my $gene_name = $fields[2];
	my $go_id = $fields[4];
	
	# top level GO categories (molecular_function, cellular_component, biological_process) - we don't need these
	next if ($go_id eq 'GO:0003674' or $go_id eq 'GO:0005575' or $go_id eq 'GO:0008150');
	
	if(exists $go_categories{$go_id}){

		my @current_array = @{$go_categories{$go_id}};
		# assuming that the order of the genes is not mixed up
		if ($current_array[-1] eq $gene_name){
			next;
		}
		else{
			push(@{$go_categories{$go_id}}, $gene_name);
		}
	}
	
	# If the GO category doesn't exist in the hash, create it. We need to look up the info from the obo file.
	else{
		
		my @go_category_info = look_up_id($go_id);
		
		my $go_term = $go_category_info[0];
		my $type = $go_category_info[1];
		
		chomp $go_term;
		# perl gets upset with the whitespace
		$go_term =~ s/ /_whitespace_/g;

		my $short_type;
		
		if($type eq 'molecular_function'){
				$short_type = 'GOMF';
		}
		elsif($type eq 'biological_process'){
				$short_type = 'GOBP';
		}
		if($type eq 'cellular_component'){
				$short_type = 'GOCC';
		}
		
		my @info = (uc($go_term), $short_type, $go_id);

		$go_categories{$go_id} = \@info;
		push(@{$go_categories{$go_id}}, $gene_name);
	}	
}
close IN;


# print out all the info and genes for each GO id
foreach my $go_id (keys %go_categories){

	my @go_info = @{$go_categories{$go_id}};
	my $term = lc(shift @go_info);
	$term =~ s/_whitespace_/ /g;
	my $identifier = join('%', uc($term), shift @go_info, shift @go_info);
	my $genes = join("\t", @go_info);

	print OUT"$identifier\t$term\t$genes\n";
}


sub look_up_id{

	my ($id) = @_;

	if(exists $go_id_description{$id}){
		
		my @info =  @{$go_id_description{$id}};
		return @info;
	}
}

# load the obo info into a hash - %go_id_description - so we can get the info such as biological_process and the name of the category.
# This is also available as a separate script - get_ontology_info_from_obo_file.pl
# We also have a separate hash to store the parents and children

sub load_obo_file_into_hash{

	my ($file) = @_;
	
	# the go id
	my $id;
	# the go category name
	my $category_name;
	# keep track of number of terms loaded 
	my $counter =0;
	
	open (IN_OBO, $file) or die "Cannot open file $file\n";
	 
	OBO_LINE: while (my $line = <IN_OBO>){

		chomp $line;
		if ($line eq '[Term]'){
			
			$counter++;
			
			$id = <IN_OBO>;
			chomp $id;
			$id =~ s/id: //;
			
			my $category_name = <IN_OBO>;
			chomp $category_name;
			$category_name =~ s/name: //;
			
			my $type = <IN_OBO>;
			chomp $type;
			$type =~ s/namespace: //;

			my @merged_name_descr = ($category_name, $type);						
			# for each id, add it to the hash
			$go_id_description{$id} = \@merged_name_descr;
			
			# there are some alternative ids, we're just going to add these in for now so there might be some duplication
			my $next_line = <IN_OBO>;
			while($next_line =~ /^alt_id/){
				$id = $next_line;
				chomp $id;
				$id =~ s/alt_id: //;
				
				$go_id_description{$id} = \@merged_name_descr;
				$next_line = <IN_OBO>;
			}	
								
			next if($line =~ /^!/);			
		}
	}
	print "\n$counter ontology categories processed from $obo_file \n";
	close IN_OBO;	
}

