#!/usr/bin/perl
use warnings;
use strict;
use CGI;
use FindBin qw($RealBin);
use File::Glob;
use HTML::Template;

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

sub show_job {

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

