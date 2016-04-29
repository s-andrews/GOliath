#!/usr/bin/perl

# This scripts connects to the ensembl api and gets all the lowest level GO terms for each gene.
# It does not iterate up and get all the parent terms (the code for this is commented out at the end of the script)  - that takes too long to do via the API so that is done in a separate script. 
# This script is designed to precede the scripts process_go_child_terms.pl and get_GO_parents.pl or get_GO_parents_and_relatives.pl.

# It requires bioperl and ensemblapi modules
#module load bioperl 
#module load ensemblapi

# USAGE: perl get_GO_child_terms_ensembl_api.pl Homo_sapiens > output_file.txt


use warnings;
use strict;
use Bio::EnsEMBL::Registry;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Location::Split;
use Bio::Location::Simple;


$|++;

system("clear") == 0 or warn "Couldn't clear screen";

my $registry = load_registry();

my $GO_adapter =   $registry->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );

die "Couldn't get GO adaptor" unless ($GO_adapter);

my $species = $ARGV[0];

my $db_adapter;

if ($species) {
  my $db_adapter = $registry->get_DBAdaptor($species,'Core');
  process_genome($db_adapter);
}
else {
  die "You need to specify a species\n";
}


sub load_registry {

  warn "Loading Registry information - please be patient\n";

  my $registry = 'Bio::EnsEMBL::Registry';

  $registry->load_registry_from_db(
				 -host => 'ensembldb.ensembl.org',
				 -user => 'anonymous'
				);

  return $registry;

}



sub process_genome {

  my $db_adapter = shift;

#  die "Adapter is $db_adapter\n";

  warn "Processing genome ".$db_adapter->species()."\n";

  my $simple_feature_adapter = $db_adapter->get_SimpleFeatureAdaptor();

  my $assembly = $db_adapter->get_adaptor('coordsystem')->fetch_all->[0]->version();
  my $species = $db_adapter->species();

  my $readable_species = $species;
  $readable_species =~ s/_/ /g;
  $readable_species =~ s/^(\w)/uc $1/e;

  # In some Ensembl species there are duplicated regions (eg the PAR in chrs X/Y).
  # If you want to see the full lengths of chromosomes then you need to pass in
  # the 4 argument form of fetch_all.  Specifying just 'chromosome' gets us separate
  # regions for the non-duplicated parts of duplicated chromosomes.

  my @chr_slices = @{$db_adapter -> get_adaptor('slice') -> fetch_all('chromosome',undef,0,1)};

  # Destroy the chr array as we're iterating since
  # the API does lazy loading which will make these
  # objects balloon in size as we process them.

  while (@chr_slices) {

    my $chr_slice = shift @chr_slices;

    process_chromosome($chr_slice,$readable_species,$assembly,$simple_feature_adapter);
  }

}


sub process_chromosome {

  my $chr_slice = shift;
  my $species = shift;
  my $assembly = shift;
  my $simple_feature_adapter = shift;

  if (length($chr_slice->seq_region_name()) > 5) {
    warn "Skipping odd looking chromsome ".$chr_slice->seq_region_name()."\n";
    return;
  }

  warn "Processing chromosome ".$chr_slice->seq_region_name()."(".$chr_slice->length()." bp) \n";

  ###### Testing only ######
  # return unless ($chr_slice->seq_region_name() eq 'Y');
  ##########################


  my @genes = @{$chr_slice -> get_all_Genes()};

  while (@genes) {
    my $gene = shift @genes;
	print $gene;
    process_gene($gene);
  }


}

sub process_gene {
  my $gene = shift;

  unless ($gene->external_name()) {
    $gene->external_name($gene->stable_id());
  }

#  warn "Found gene ".$gene->external_name()." of type ".$gene->biotype()."\n";

  my $genetype = 'gene';

  if ($gene->biotype() eq 'pseudogene') {
    $genetype = 'pseudogene';
  }

  my @xrefs = @{$gene->get_all_DBEntries()};

  my @transcripts = @{$gene -> get_all_Transcripts()};

  while (@transcripts) {
    my $transcript = shift @transcripts;

    if ($gene -> biotype() !~ /pseudogene/) {
      process_transcript($transcript,$gene);
    }
  }

}

sub process_transcript {

  my $transcript = shift;
  my $gene = shift;

  unless ($transcript->external_name()) {
    $transcript->external_name($gene->external_name());
  }

#  warn "Found transcript ".$transcript->external_name()." of type ".$transcript->biotype()."\n";

  process_translation($transcript,$gene);
}


sub process_translation {

  my $transcript = shift;
  my $gene = shift;
  my $translation = $transcript->translation();

  # Now see if we have a translation
  my @translateable_exons = @{$transcript->get_all_translateable_Exons()};
  my $location;

  if (@translateable_exons == 0) {
    # No translation
    return;
  }


#  warn "Found translation ".$transcript->external_name()."\n";

  my @xrefs = @{$translation->get_all_DBEntries()};

  while (@xrefs) {
    my $xref = shift @xrefs;
    next unless($xref); # Shouldn't ever be empty, but we've seen this happen.
    next unless ($xref->dbname());

    if ($xref -> dbname() eq 'GO') {
      my $term = $GO_adapter->fetch_by_accession($xref->display_id());

	  next unless ($term); # Shouldn't happen, but can.
#      $feature -> add_tag_value(db_xref => $xref->display_id()." ".$term->name." (".$term->namespace().")");

      print_terms($gene,$term);


    }

#    warn "Xref of ".$xref->display_id()." in ".$xref->dbname()."\n";
  }

}

sub print_terms {
  my ($gene,@terms) = @_;

  foreach my $term(@terms) {
  
	# to print out gene symbols
    #print join("\t",($gene->external_name(),$term->accession(),$term->namespace(),$term->name())),"\n";
	# to print out ensembl ids (display_id)
	print join("\t",($gene->display_id(),$term->accession(),$term->namespace(),$term->name())),"\n";
  
 # uncomment the following if we want to get all the parent terms - it takes a long time.
  # my $parents = $term->parents();
   # if (@$parents) {
   #   print_terms($gene,@$parents);
   # }

  }
}
