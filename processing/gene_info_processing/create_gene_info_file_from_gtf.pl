#!/usr/bin/perl

# This script parses a gtf file and uses it to look up information such as GC content from the appropriate genome

# download gtf from ensembl:
#     ftp://ftp.ensembl.org/pub/current_gtf/

# /bi/apps/clusterflow/_devel/genomes.config

use warnings;
use strict;
$|++;
use Getopt::Long;
use Cwd;
use feature 'say';

my $gtf_file;
my $genome;
my $genome_folder;
my $exclude_biotypes;
my $include_biotypes;
my $list_biotypes;
#my @biotypes_to_exclude;
my %biotype_exclusion_count;
my %biotype_inclusion_count;
my %gene_biotypes;
my %biotype_families;
my $biotype_family_file;
my $help;
my $outfile_suffix = '_gene_info.txt';

my $config_result = GetOptions(
"gtf=s" => \$gtf_file,
"genome=s"  => \$genome,
"genome_folder=s"  => \$genome_folder,
"list_biotypes" => \$list_biotypes,
"include_biotypes=s" => \$include_biotypes,
"exclude_biotypes=s" => \$exclude_biotypes,
"biotype_family_file=s" => \$biotype_family_file,
"help" => \$help
);
die "Could not parse options" unless ($config_result);


if ($help) {
    print while (<DATA>);
    exit;
}

unless((defined($gtf_file)) && ($gtf_file =~ /.gtf$|.gtf.gz$/)){
	die "\ngtf file must be supplied using --gtf file_path\n";
}	

unless(defined($biotype_family_file)){
	$biotype_family_file = "biotype_families.txt";
}	

read_in_biotype_families();


#Check if the read file is compressed and open accordingly
if($gtf_file =~ /\.gz$/){
	print "\nusing gtf file: $gtf_file\n";
	print "----------------------------------\n";
	open (GTF_IN,"zcat \'$gtf_file\' |") or die "Can't read $gtf_file: $!";
}
else{
	print "\nusing gtf file: $gtf_file\n";
	#print "----------------------------------\n";
	open (GTF_IN,$gtf_file) or die "Can't read $gtf_file: $!";
}

if($list_biotypes){

	print "\ngene biotypes in gtf file...\n";
	print "------------------------------\n";
	
	while (my $line = <GTF_IN>){
		
		# skip the header lines
		if($line =~ /^#/){
			next;
		}	
		
		my @line_info = split(/\t/, $line);
		
		if($line_info[2] eq "gene"){ 
		
			my @gene_info = split(/;/, $line_info[8]);
			
			my $biotype = $gene_info[4];
			$biotype =~ s/gene_biotype "//;
			$biotype =~ s/\"//g;
			$biotype =~ s/ //g;
			
			if (exists $gene_biotypes{$biotype}){
				$gene_biotypes{$biotype} ++;
			}	
			else{
				$gene_biotypes{$biotype} = 1;
			}
		}	
	}
	
	# print the ordered hash
	foreach my $gene_biotype(sort {$gene_biotypes{$b} <=> $gene_biotypes{$a}} keys %gene_biotypes){

		print $gene_biotype, "\t", $gene_biotypes{$gene_biotype},"\n";
	}	
	exit;
}	


if(defined($exclude_biotypes)){
	my @biotypes_to_exclude = split(" ", $exclude_biotypes);
	chomp @biotypes_to_exclude;
	
	$outfile_suffix = '_excluding_'.$exclude_biotypes.'_biotypes.txt';
	$outfile_suffix =~ s/ /_/g;

	foreach my $biotype_to_exclude(@biotypes_to_exclude){
		if (exists $biotype_exclusion_count{$biotype_to_exclude}){
			warn "skipping biotype $biotype_to_exclude as it appears to be duplicated\n";
		}
		else{
			$biotype_exclusion_count{$biotype_to_exclude} = 1;
			print "excluding gene biotype $biotype_to_exclude\n";
		}	
	}
	print "----------------------------------\n";
}

elsif(defined($include_biotypes)){
	my @biotypes_to_include = split(" ", $include_biotypes);
	chomp @biotypes_to_include;
	
	$outfile_suffix = "_only_including_".$include_biotypes.'_biotypes.txt';
	$outfile_suffix =~ s/ /_/g;

	foreach my $biotype_to_include(@biotypes_to_include){
		if (exists $biotype_inclusion_count{$biotype_to_include}){
			warn "skipping biotype $biotype_to_include as it appears to be duplicated\n";
		}
		else{
			$biotype_inclusion_count{$biotype_to_include} = 1;
			print "including gene biotype $biotype_to_include\n";
		}	
	}
	print "----------------------------------\n";	
}

my %chromosomes;
my $parent_dir = getcwd;

# check whether a genome folder has been supplied
if(defined($genome_folder)){
	print "\nusing genome_folder: $genome_folder\n";
	read_genome_into_memory();
}
elsif(defined($genome)){
	get_genome_path();
	read_genome_into_memory();
}	
else{
	warn "\nGenome folder not supplied, info file will still be created but will not contain GC content.\nTo supply a genome location use the option --genome\n";
}

my $outfile = $gtf_file;
$outfile =~ s/.gtf$//;
$outfile =~ s/.gtf.gz$//;
$outfile = $outfile.$outfile_suffix;

open(OUT, '>', $outfile) or die;

my $gene_counter = 0;
my $line_counter = 0;

my $gene_id;
my $gene_name;
my $chr;
my $start;
my $end;
my $strand;
my $biotype;
my $biotype_family;
my $length;
my $GC_content;
my $no_of_transcripts;
my $sequence;


# print the header line
my $header = "gene_id\tgene_name\tchromosome\tstart\tend\tstrand\tbiotype\tbiotype_family\tlength\tGC_content\tno_of_transcripts";
print OUT "$header\n"; 
 
 
print "\nprocessing gtf file...\n";

LINE_LOOP: while (my $line = <GTF_IN>){

	$line_counter++;

	# skip the header lines
	if($line =~ /^#/){
		next;
	}	
	
	#print "whole line: $line\n";
	
	my @line_info = split(/\t/, $line);
	
	# if it's not a proper chromosome in a supplied genome assembly then we don't want to use it
    # but if a genome assembly hasn't been provided we'll just work with what is in the gtf file, even if they are pseudo chromosomes...
	if(defined($genome_folder)){
	
		unless(exists($chromosomes{$line_info[0]})){
			#print("chromosome not found: $line_info[0] \n");
			next LINE_LOOP;
		}	
	}	
	
	if($line_info[2] eq "gene"){
	
		#remove any whitespace at end
		$line_info[8] =~ s/\s+$//;
	
		my @gene_info = split(/;/, $line_info[8]);
		
		my %attributes;
		
		# the attributes in the gtf file are not all the same in the different species - some fields are missing
		# load them into the hash here so we can check the biotypes, then check them all after the printing to file has been done.
		foreach my $attribute_info(@gene_info){
			
			#remove any whitespace at start
			$attribute_info =~ s/^\s+//;
			
			#print("\nattribute info: = $attribute_info");
			my @attribute = split(/\s+/, $attribute_info);
			#my $attribute_name = $attribute[0];
			#print("\nattribute name = $attribute_name");
			my $attribute_value = $attribute[1];
			
			# strip the quotes and any whitespace
			$attribute_value =~ s/\"|\s+//g;
			#$attribute_value =~ s/\s+//g;
			#print("\nattribute value = $attribute_value");
						
			$attributes{$attribute[0]} = $attribute_value;
			
		}
		
	#	my $biotype_to_check =  $gene_info[4];
	#	$biotype_to_check =~ s/gene_biotype "|\"| //g;
		my $biotype_to_check;
		if(exists $attributes{gene_biotype}){
			$biotype_to_check = $attributes{gene_biotype};
		}
		else{
			$biotype_to_check = "";
		}
		#print("\nbiotype to check = $biotype_to_check");
				
		# work out what I'm doing with the exclusion counts - where we actually want to do the counting
		
		if(defined($exclude_biotypes)){
			if (exists $biotype_exclusion_count{$biotype_to_check}){
				$biotype_exclusion_count{$biotype_to_check} ++;
				next LINE_LOOP;
			}	
		}	
		
		if(defined($include_biotypes)){
			unless(exists $biotype_inclusion_count{$biotype_to_check}){	
				next LINE_LOOP;
			}
		}		
	
		# write out the gene info unless it's the first one
		if($gene_counter > 0){
		
			#unless($no_of_transcripts > 0){
			#	warn "Couldn't find any transcripts for $gene_name\n";
			#}	
		
			my $output_line = join("\t", $gene_id, $gene_name, $chr, $start, $end, $strand, $biotype, $biotype_family, $length, $GC_content, $no_of_transcripts);
			print OUT "$output_line\n";
			
			if($gene_counter % 2000 == 0){
				print("\n$gene_counter genes processed...");
			}
			# reset no of transcripts
			$no_of_transcripts = 0;
		}	
	
		# get the information that we want
		$biotype = $biotype_to_check;
		
		if(exists $attributes{gene_id}){
			$gene_id = $attributes{gene_id};
		}
		else{
			$gene_id = "";
		} 
	
		if(exists $attributes{gene_name}){
			$gene_name = $attributes{gene_name};
		}
		else{
			$gene_name = "";
		} 
	
		# these should always be included in the file
		$chr = $line_info[0];
		$start = $line_info[3];
		$end = $line_info[4];
		$strand = $line_info[6];
		$length = $end-$start;
		
		if(defined($genome_folder) && length($chromosomes{$chr}) > ($start+$length)){
			$sequence = substr($chromosomes{$chr},$start,$length);
			#print("\n sequence: $sequence");
			$GC_content = get_GC_content();
			#print("\n GC content: $GC_content");
		}
		
		else{
			$GC_content = "";
		}
		
		if(exists($biotype_families{$biotype})){
			$biotype_family = $biotype_families{$biotype};
		}
		else{
			$biotype_family="";
		}	
		
		$gene_counter++;
	}
	elsif($line_info[2] eq "transcript"){	
		
		my @transcript_gene_info = split(/;/, $line_info[8]);
		my $transcript_gene_id = $transcript_gene_info[0];
		$transcript_gene_id =~ s/gene_id "//;
		$transcript_gene_id =~ s/\"//g;
		
		# check whether the transcript gene id is the same as the current gene id - if not the gene has been skipped due to biotype
		unless(defined($gene_id)){
			next LINE_LOOP;
		}	
		unless($transcript_gene_id eq $gene_id){
			next LINE_LOOP;
		}	
		$no_of_transcripts++;
	}	
}	

# write out last one, check it's not empty first.
if(defined($gene_id)){
	my $output_line = join("\t", $gene_id, $gene_name, $chr, $start, $end, $strand, $biotype, $biotype_family, $length, $GC_content, $no_of_transcripts);
	print OUT "$output_line\n";
	$gene_counter++;
}
print "\nCreated output file named $outfile containing $gene_counter genes.\n"; 

if(defined($exclude_biotypes)){

	print "\nexcluded  the following number of genes based on biotype\n";
	foreach my $gene_biotype(keys %biotype_exclusion_count){

		print $gene_biotype, " ", $biotype_exclusion_count{$gene_biotype},"\n";
	}
}

sub get_GC_content{

	my $Gcount = 0;
	while ($sequence =~ /G/g) {
		$Gcount++;
	}

	my $Ccount = 0;
	while ($sequence =~ /C/g) {
		$Ccount++;
	}

	my $seq_length = length($sequence);
	my $GC_proportion = (($Gcount + $Ccount)/$seq_length);
	my $GC_rounded = sprintf("%.3f", $GC_proportion);
	return $GC_rounded;
}

sub get_genome_path{

#	my $genome = "GRCm38";
	my $cf_genome_config_file = "/bi/apps/clusterflow/0.3_devel/genomes.config";
	open (GENOME_IN,$cf_genome_config_file) or die "Can't read $cf_genome_config_file: $!";

	while (my $line = <GENOME_IN>){
			
		# only intereseted in the genome path lines
		if($line =~ /^\@genome_path/){

			my @line_info = split(/\t/, $line);
			
			if($line_info[1] eq $genome){
			
				$genome_folder = $line_info[2];
				print("\nUsing genome from $genome_folder\n");
				return;
			}
		}	
	}	
}

sub read_genome_into_memory{

  ## reading in and storing the specified genome in the %chromosomes hash
  chdir ($genome_folder) or die "Can't move to $genome_folder: $!";
  warn "Now reading in and storing sequence information of the genome specified in: $genome_folder\n\n";

  my @chromosome_filenames =  <*.fa>;

  ### if there aren't any genomic files with the extension .fa we will look for files with the extension .fasta
  unless (@chromosome_filenames){
    @chromosome_filenames =  <*.fasta>;
  }
  unless (@chromosome_filenames){
    die "The specified genome folder $genome_folder does not contain any sequence files in FastA format (with .fa or .fasta file extensions)\n";
  }

  foreach my $chromosome_filename (@chromosome_filenames){

    # skipping the tophat entire mouse genome fasta file
    # next if ($chromosome_filename eq 'Mus_musculus.NCBIM37.fa');

    open (CHR_IN,$chromosome_filename) or die "Failed to read from sequence file $chromosome_filename $!\n";
    ### first line needs to be a fasta header
    my $first_line = <CHR_IN>;
    chomp $first_line;
    $first_line =~ s/\r//; # removing /r carriage returns

    ### Extracting chromosome name from the FastA header
    my $chromosome_name = extract_chromosome_name($first_line);
	
    my $sequence;
    while (<CHR_IN>){
      chomp;
      $_ =~ s/\r//; # removing /r carriage returns

      if ($_ =~ /^>/){
	### storing the previous chromosome in the %chromosomes hash, only relevant for Multi-Fasta-Files (MFA)
	if (exists $chromosomes{$chromosome_name}){
	  warn "chr $chromosome_name (",length $sequence ," bp)\n";
	  die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name!\n";
	}
	else {
	  if (length($sequence) == 0){
	    warn "Chromosome $chromosome_name in the multi-fasta file $chromosome_filename did not contain any sequence information!\n";
	  }
	  warn "chr $chromosome_name (",length $sequence ," bp)\n";
	  $chromosomes{$chromosome_name} = $sequence;
	}
	### resetting the sequence variable
	$sequence = '';
	### setting new chromosome name
	$chromosome_name = extract_chromosome_name($_);
      }
      else{
	$sequence .= uc$_;
      }
    }

    if (exists $chromosomes{$chromosome_name}){
      warn "chr $chromosome_name (",length $sequence ," bp)\t";
      die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name.\n";
    }
    else{
      if (length($sequence) == 0){
	warn "Chromosome $chromosome_name in the file $chromosome_filename did not contain any sequence information!\n";
      }
      warn "chr $chromosome_name (",length $sequence ," bp)\n";
      $chromosomes{$chromosome_name} = $sequence;
    }
  }
  warn "\n";
  chdir $parent_dir or die "Failed to move to directory $parent_dir\n";
} # end read_genome_into_memory

sub extract_chromosome_name {
  ## Bowtie extracts the first string after the initiation > in the FASTA file, so we are doing this as well
  my $fasta_header = shift;
  if ($fasta_header =~ s/^>//){
    my ($chromosome_name) = split (/\s+/,$fasta_header);
    return $chromosome_name;
  } # end if
  else{
    die "The specified chromosome ($fasta_header) file doesn't seem to be in FASTA format as required!\n";
  }
} # end extract_chromosome_name

sub read_in_biotype_families{
	
	warn "\nNow reading in and storing gene biotype information specified in: $biotype_family_file\n\n";
	# check whether file exists......
	unless(-e $biotype_family_file){
		warn "No biotype family file found, the processing can be run without it, but a file can easily be created by entering the mysql command in the help section\n";
		#warn ('mysql -uanonymous -P3306 -hensembldb.ensembl.org -Densembl_production_73 -e "select distinct(name),biotype_group from biotype where db_type like '%core%' and is_current=1 order by biotype_group,name;"  > biotype_families.txt')
		return;
	}	
	open (BIOTYPE_IN,$biotype_family_file) or die "Failed to read from biotype file $biotype_family_file $!\n";		
	# to check the families
	my %biotype_family_count;
	
	my @header = split(/\t/, <BIOTYPE_IN>);
	chomp @header;
	
	unless(($header[0] eq "name") && ($header[1] eq "biotype_group")){
		warn "first line of biotype family file is not as expected\n";
	}	
		
	
	while (<BIOTYPE_IN>){
		
		
		chomp $_;
		
		my @line = split(/\t/, $_);
		
		if (exists $biotype_families{$line[0]}){
		
			warn "$line[0] already ezists\n";
		}	
		else{							
			$biotype_families{$line[0]} = $line[1];
			
			# to check the biotype family file looks ok
			if (exists $biotype_family_count{$line[1]}){
			
				$biotype_family_count{$line[1]}++;
			}
			else{
				$biotype_family_count{$line[1]} = 0;
			}	
		}
	}
	print "=====================================\n";
	print "biotype families read in from file\n";
	print "-------------------------------------\n";
	
	foreach my $biotype_family(keys %biotype_family_count){
		print "$biotype_family\t$biotype_family_count{$biotype_family}\n";
	}
	print "=====================================\n";
	close BIOTYPE_IN;	
}


__DATA__

========================================================================================================================

 Perl script that parses a gtf file and uses it to look up information such as GC content from the appropriate genome.

 Usage: create_gene_info_from_file.pl --gtf file.gtf --genome GRCm38 --exclude_biotypes "x y z"

========================================================================================================================

Options:

  --gtf                   gtf or gtf.gz file containing gene and transcript information (required)

  --genome                genome to use e.g. GRCm38 (optional, GC content column will be blank if this isn't specified)
                          uses the clusterflow config file to get the genome paths
						  run cf --list_genomes to see available genomes 

  --genome_folder         path to genome folder - another way of accessing the genome. This isn't needed if --genome is specified.

  --biotype_family_file   path to file containing biotypes and family they belong to
                          see command at the end of the help section

  --exclude_biotypes      gene biotypes to exclude - string of biotypes enclosed within quotes and separated by whitespace
                            e.g. --exclude_biotypes "processed_pseudogene another_type protein_coding"

  --include_biotypes      gene biotypes to include - only these biotypes will be included. 
                          use either this option or the exclude biotypes options, do not use both in the same command
                          format is the same as for exclude_biotypes

  --list_biotypes         lists all the gene biotypes along with the number of genes and exits

  --help                  display help and exit

========================================================================================================================

output file should contain 
gene_id  gene_name  chromosome  start  end  strand  biotype  length  GC_content  no_of_transcripts

gtf files can be downloaded from ftp://ftp.ensembl.org/pub/current_gtf/....
This script assumes the files are in this format. 


to get the biotype families - this isn't species specific - this exact command will work
mysql -uanonymous -P3306 -hensembldb.ensembl.org -Densembl_production_73 -e "select distinct(name),biotype_group from biotype where db_type like '%core%' and is_current=1 order by biotype_group,name;"  > biotype_families.txt

