#!/usr/bin/perl

# This takes the processed child terms and finds all the parent categories. To be used after the process_go_child_terms.pl script. 
# This is just looking up the parents i.e. is_a relationships.

# usage perl get_GO_parents.pl xxxx_go_terms.txt go-basic.obo

use warnings;
use strict;
use Data::Dumper;

my $gmt_kids_file = $ARGV[0];
my $obo_file = $ARGV[1];
#my $obo_file = "go-basic.obo";

# key is child, values are parents
my %child_parents;

# used to store all the go ids and descriptions.
my %go_id_description;

my $id; 
my @parent_ids;
my $counter = 0;
my $orphan_count;

# process the obo file
open (IN, $obo_file) or die "Cannot open file $obo_file\n";

OBO_LINE: while (my $line = <IN>){

	chomp $line;
	if ($line eq '[Term]'){
		
		$counter++;
		
		if(defined($id)){

			if(scalar(@parent_ids) == 0){
			$orphan_count++;
				#warn "no parent found for $id \n";				
			}	
		
			else {
				
				my @parent_ids_to_add = @parent_ids;
				$child_parents{$id} = \@parent_ids_to_add;
				
				# clear the parent ids
				@parent_ids = ();
			}	
		}
		
		$id = <IN>;
		chomp $id;
		$id =~ s/id: //;
		my $category_name = <IN>;
		chomp $category_name;
		$category_name =~ s/name: //;
		my $type = <IN>;
		chomp $type;
		$type =~ s/namespace: //;
		
		my $processed_type = process_type($type);		
			
		my $merged_name_descr = join("%", $category_name, $processed_type);
		
		# for each id, add it to the hash
		$go_id_description{$id} = $merged_name_descr;
		
		#print "\nnext id is $id \n";
	}
	elsif($line =~ m/^is_a:/){
	
		# This is assuming there is always whitespace after is_a: . I checked this (Nov 14) and it was always the case.
		$line =~ s/^is_a: //;
		my @parent_info = split(/ ! /, $line);
		my $parent_id = $parent_info[0];
		
		# In the go-basic.obo there are a couple of lines that have is_a: regulates ! regulates
		if ($parent_id =~ m/^GO:/){
			push (@parent_ids, $parent_id);
		}
		else {
			warn "Did not recognise $parent_id as a GO id \n";
		}	
	}
}

#print Dumper %child_parents;
print "\n============================================================\n";
print "\n$counter ontology categories processed from $obo_file \n";
print "\nno parents found for $orphan_count categories \n";
print "\n============================================================\n";
close IN;



open (IN, $gmt_kids_file) or die "Cannot open file $gmt_kids_file\n";

my $outfile = $gmt_kids_file;
$outfile =~ s/.txt//;
$outfile =~ s/just_children.//;
$outfile = "$outfile.including_parents.gmt";
open (OUT, ">", $outfile) or die  "Couldn't open output file: $! \n";
#open (OUT_LOG, ">", "log.txt") or die  "Couldn't open output file: $! \n";

my %gmt_file;
my $child_counter;
my $parent_counter;
my $line_counter;
my $go_category_counter;

$|++;

# read in the whole gmt file 
while (my $line = <IN>){

	$line_counter++;
	chomp $line;
	my @fields = split(/\t/, $line);
	my @info = split(/%/, $fields[0]);
	my $go_id = $info[2];
	chomp $go_id;
	
	$gmt_file{$go_id} = \@fields;
	
}

print "loaded $line_counter categories\n";
$line_counter = 0;

foreach my $go_id (keys %gmt_file){

	$line_counter++;
	my @gmt_line = @{$gmt_file{$go_id}};
	
	my @info = split(/%/, shift @gmt_line);

	if($go_id ne $info[2]){
		warn "something's gone wrong, the GO IDs don't match\n";
	}
	else{	
		# remove the description
		my $descr = shift @gmt_line;
		my @genes = @gmt_line;
		chomp(@genes);

		# now go through each category and populate the parents
		my @parents = get_parents($go_id, %child_parents);
		my $count = scalar(@parents);
		
		#print "found $count parents for $go_id\n";
		
		while(@parents){

			my @all_parents = ();
		
			foreach my $parent(@parents){
				
				#print "parent id = $parent\n";
				
				# create the go category if it doesn't already exist - look up the information for the GO id in the %go_id_description hash
				if(exists $gmt_file{$parent} != 1)  {

					if (exists $go_id_description{$parent}){
									
						my ($go_term, $processed_type) = get_ontology_info_using_id($parent, %go_id_description);						
						
						%gmt_file = add_new_go_category($processed_type, $go_term, $parent, %gmt_file);
						
						#print "adding new category for id $parent and term $go_term\n";
						
						$go_category_counter++;
					}
					else {
						warn "\nCouldn't find parent id $parent in the obo file \n\n";
					}	
				}
				# check that it does now exist and add the genes - we don't want the genes to be duplicated hence the check for whether they alrady exist.
				if(exists $gmt_file{$parent}){
					
					#print "It exists so we're trying to add the genes to parent $parent\n";
					
					my @original_gmt_line = @{$gmt_file{$parent}};
					
					my %hash;
					$hash{$_}++ for (@original_gmt_line);
					
					foreach my $gene(@genes){
					
						next if (exists $hash{$gene});
						
						push(@{$gmt_file{$parent}}, $gene);
					}
										
					$parent_counter++;				
				}
			
				# get the next parents	
				my @next_generation = get_parents($parent, %child_parents);				
				push(@all_parents,  @next_generation);
							
			}	
			# deduplicate the parent ids
			my %hash   = map { $_, 1 } @all_parents;
			my @unique = keys %hash;
			@parents = @unique;
			#@parents = @all_parents;

			#print "\n\nchecking for duplication\n";
			#foreach my $test(@parents){
			#	print "parent found = $test\n";
			#}				
		}
	}
	if(($line_counter % 100) == 0){
		print "processed $line_counter child terms, added $go_category_counter parent categories\n";
	}	
}	

warn "completed processing, now attempting to write out results to $outfile \n";

foreach my $go_id (keys %gmt_file){

	my @go_info = @{$gmt_file{$go_id}};
	my $identifier = shift @go_info;
	my $term = shift @go_info;
	my $genes = join("\t", @go_info);

	print OUT "$identifier\t$term\t$genes\n";
}



sub get_ontology_info_using_id{

	my ($go_id, %go_id_description) = @_;
	
	my $ontology_info = $go_id_description{$go_id};
	
	my @ontology_info_fields = split(/%/, $ontology_info);	

	my $go_term = $ontology_info_fields[0];
	
	my $processed_type = $ontology_info_fields[1];
	
	return($go_term, $processed_type);
}

	
sub add_new_go_category{
	
	my ($processed_type, $go_term, $go_id, %go_categories) = @_;
	
	# we want the information about the ontology category to be the first items in the array, then we add the genes afterwards.
	my $identifier = join('%', uc($go_term),  $processed_type, $go_id);
	
	my @info = ($identifier, $go_term);

	$go_categories{$go_id} = \@info;
	
	return %go_categories;
}	
	
sub get_parents {

	my ($go_id, %child_parent_hash) = @_;
	
	#print "child id in subroutine = $go_id\n";
	
	my @filtered_parents = ();
	
	if(exists $child_parent_hash{$go_id}){
	
		my @parents =  @{$child_parent_hash{$go_id}};
		
		foreach my $parent_id(@parents){
	
			#print "parent id in subroutine = $parent_id\n";
	
			next if ($parent_id eq 'GO:0003674' or $parent_id eq 'GO:0005575' or $parent_id eq 'GO:0008150');
			next unless ($parent_id =~ /^GO:/);
			
			#print "made it through\n";
			
			push(@filtered_parents, $parent_id);
		}
	}
	return @filtered_parents;
}	



sub process_type{

	my ($category_type) = @_;
	my $string;
	
	if($category_type eq 'molecular_function'){
		$string = 'GOMF';
		return $string;
	}
	elsif($category_type eq 'biological_process'){
		$string = 'GOBP';
		return $string;
	}
	elsif($category_type eq 'cellular_component'){
		$string = 'GOCC';
		return $string;
	} 
	else{
		warn "category type could not be identified for $category_type\n";
	}	
}	




	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


