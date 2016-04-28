#!/usr/bin/perl
use warnings;
use strict;
use CGI;
use FindBin qw($RealBin);
use File::Glob;
use HTML::Template;
use CGI::Carp qw(fatalsToBrowser);

#######################################################################
# Copyright Simon Andrews (simon.andrews@babraham.ac.uk) 2016
#
# This file is part of GOliath.
#
# GOliath is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# GOliath is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GOliath.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################

my $q = CGI -> new();

my $job_id = $q -> param("job_id");

if ($job_id) {
    show_job($job_id);
}
elsif ($q -> param("submit")) {
    process_submission();
}
else {
    show_home();
}

sub show_home {

    my $template = HTML::Template -> new(filename => "$RealBin/../templates/home.html");
    
    my @species = list_species();

    unless (@species) {
	print_bug("Couldn't find any valid godata");
    }
    
    my @species_templates;

    foreach my $species (@species) {
	push @species_templates, {SPECIES => $species};
    }

    $template -> param(SPECIES => \@species_templates);

    print $template -> output();


}

sub print_bug {
    my ($message) = @_;

    my $template = HTML::Template -> new(filename => "$RealBin/../templates/bug.html");
    $template -> param(MESSAGE => $message);
    print $template->output();
    die $message;
}

sub print_error {

    my ($message) = @_;

    my $template = HTML::Template -> new(filename => "$RealBin/../templates/error.html");
    $template -> param(MESSAGE => $message);
    print $template->output();
    warn $message;
    exit;

}

sub process_submission {

    # We need to collect the information from the form and check that it works.
    my $species = $q -> param("species");
    my $list_type = $q -> param("list_type");
    my $gene_list_text = $q -> param("gene_list");
    my $background_list_text = $q -> param("background_list");

    unless ($species) {
	print_bug("No species when submitting form");
    }
    my @species = list_species();
    my $valid_species;

    foreach my $known_species (@species) {
	if ($known_species eq $species) {
	    $valid_species = $known_species;
	    last;
	}
    }

    unless ($valid_species) {
	print_bug("Couldn't find species '$species' on the server");
    }

    unless ($list_type eq 'Ordered' or $list_type eq 'Unordered') {
	print_bug("Uknown list type '$list_type'");
    }

    my @gene_list_genes = get_genes($gene_list_text,$valid_species);
    my @background_list_genes = get_genes($background_list_text,$valid_species);

    unless (@gene_list_genes) {
	print_bug("Found no gene list genes");
    }

    if ($list_type eq "Ordered") {
	@background_list_genes = ();
    }

    my $job_id = generate_job_id();

    print $q->redirect("goliath.cgi?job_id=$job_id");

}

sub generate_job_id {

    my @letters = ('A'..'Z','a'..'z',0..9);

    while (1) {
	my $code;

	for (1..20) {
	    $code .= $letters[int rand(scalar @letters)];
	}

	if (-e "$RealBin/../../jobs/$code") {
	    warn "Code $code already exists";
	    next;
	}

	return $code;
	

    }


}



sub get_genes {
    my ($text,$species) = @_;

    # Eventually we'll be more clever about this and will validate and deduplicate
    # the lists, but for now we'll just split them and be done with it.

    $text =~ s/$\s+//g;
    $text =~ s/\s+$//g;

    my @genes = split(/\s*[\r\n]+\s*/,$text);

    return @genes;


}


sub show_job {
    my ($job_id) = @_;
    die("Show job not implemented yet for $job_id");
}


sub list_species {
    my @valid_species;

    chdir("$RealBin/../../godata") or print_bug("Couldn't move to godata folder: $!");

    my @species = File::Glob::bsd_glob("*");


    foreach my $species (@species) {
	next unless (-d $species);

	my @versions = File::Glob::bsd_glob("$species/*");

	foreach my $version (@versions) {
	    next unless (-d $version);
	    # TODO: Check for the correct files inside the directory
	    push @valid_species,$version;
	}
    }
    return @valid_species;

}

