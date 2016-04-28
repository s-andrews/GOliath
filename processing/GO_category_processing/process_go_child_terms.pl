#!/usr/bin/perl

# processes the gene ontology data downloaded using the Ensembl API (using the script go_test_only_child_process.pl)

# input file format eg.
# OR4F17  GO:0050877      biological_process      neurological system process

# produces output like the gmt files used at the Broad and Bader lab e.g.
# TRANSPORTER COMPLEX%GOCC%GO:1990351     transporter complex     ccb-2   ccb-1

# It took too long to get all the parent terms using the API so that's now done after this using the script get_GO_parents.pl or get_GO_parents_realatives.pl and using the go-basic.obo file.

use warnings;
use strict;

my $file = $ARGV[0];

open (IN, $file) or die "Cannot open file $file\n";

my $outfile = $file;
$outfile =~ s/.txt//;
$outfile = "$outfile.go_terms.txt";
open (OUT, ">", $outfile) or die  "Couldn't open output file: $! \n";
my %go_categories;

while (my $line = <IN>){

	my @fields = split(/\t/, $line);
	my $gene_name = $fields[0];
	my $go_id = $fields[1];

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
	else{

		my $type = $fields[2];
		my $go_term = $fields[3];
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

foreach my $go_id (keys %go_categories){

	my @go_info = @{$go_categories{$go_id}};
	my $term = lc(shift @go_info);
	$term =~ s/_whitespace_/ /g;
	my $identifier = join('%', uc($term), shift @go_info, shift @go_info);
	my $genes = join("\t", @go_info);

	print OUT"$identifier\t$term\t$genes\n";
}
