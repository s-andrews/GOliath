#!/usr/bin/perl
use warnings;
use strict;
use CGI;
use FindBin qw($RealBin);
use File::Glob;
use HTML::Template;
#use CGI::Carp qw(fatalsToBrowser);

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

my $config = read_config();

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

    $template -> param(VERSION => $config->{VERSION},
		       ADMIN_EMAIL => $config->{ADMIN_EMAIL},
		       ADMIN_NAME => $config->{ADMIN_NAME});


    print $template -> output();


}

sub print_bug {
    my ($message) = @_;

    my $template = HTML::Template -> new(filename => "$RealBin/../templates/bug.html");
    $template -> param(MESSAGE => $message);
    $template -> param(VERSION => $config->{VERSION},
		       ADMIN_EMAIL => $config->{ADMIN_EMAIL},
		       ADMIN_NAME => $config->{ADMIN_NAME});
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

    # Making the job id should have created a folder in the jobs
    # folder for us.

    chdir ("$config->{JOB_FOLDER}/$job_id") or print_bug("Failed to move to job folder $job_id: $!");

    # Now we can save the files
    open (OUT,'>','config.txt') or print_bug("Failed to write to config.txt: $!");
    print OUT "type\t$list_type\n";
    print OUT "species\t$RealBin/../../godata/$valid_species\n";
    close OUT or print_bug("Failed to write to config.txt: $!");

    open (OUT,'>','gene_list.txt') or print_bug("Failed to write to gene_list.txt: $!");
    foreach my $gene (@gene_list_genes) {
	print OUT $gene,"\n";
    }
    close (OUT) or print_bug("Failed to write to gene_list.txt: $!");

    open (OUT,'>','background_list.txt') or print_bug("Failed to write to background_list.txt: $!");
    foreach my $gene (@background_list_genes) {
	print OUT $gene,"\n";
    }
    close (OUT) or print_bug("Failed to write to background_list.txt: $!");

    # Now we need to launch the actual analysis process.  For now we're just going to launch each
    # one as it comes in and see if we cope.  If it gets too bad we might have to institute a
    # queueing system.

    # We don't want to wait for the child so we'll detach from it
    $SIG{CHLD} = 'IGNORE';

    my $pid = fork();

    if ($pid) {
	# We're the child and we need to start the analysis

	# First we'll write out pid into a file in the run
	# folder so that the results tracker can tell if we've
	# died.

	open (PID,'>','pid.txt') or die "Failed to write to pid file : $!";

	print PID $pid;
	close PID or die "Failed to write to pid file: $!";

#	exec("Rscript $RealBin/../../analysis/GO_analysis.r \"$config->{JOB_FOLDER}/$job_id\" > log.txt 2>errors.txt");

	system("Rscript $RealBin/../../analysis/GO_analysis.r \"$config->{JOB_FOLDER}/$job_id\" > log.txt 2>errors.txt &");
	exit(0);

    }


    print $q->redirect("goliath.cgi?job_id=$job_id");

}

sub generate_job_id {

    my @letters = ('A'..'Z','a'..'z',0..9);

    while (1) {
	my $code;

	for (1..20) {
	    $code .= $letters[int rand(scalar @letters)];
	}

	if (-e "$config->{JOB_FOLDER}/$code") {
	    warn "Code $code already exists";
	    next;
	}

	unless (mkdir("$config->{JOB_FOLDER}/$code")) {
	    # The chances of generating the same code at the same time
	    # are pretty small so we'll assume this is a bug
	    print_bug("Failed to make job folder for $code: $!");
	}

	return $code;
	

    }


}



sub get_genes {
    my ($text,$species) = @_;

    # Don't do anything if there's nothing in the text;
    return () unless ($text);

    # Eventually we'll be more clever about this and will validate and deduplicate
    # the lists, but for now we'll just split them and be done with it.

    $text =~ s/^\s+//g;
    $text =~ s/\s+$//g;

    my @genes = split(/\s*[\r\n]+\s*/,$text);

    return @genes;


}


sub show_job {
    my ($job_id) = @_;

    # If the job is complete we should find a flag file called finished.flag in
    # the run folder

    my $exists = -e "$config->{JOB_FOLDER}/$job_id";

    my $complete = 0;

    if ($exists) {
	chdir ("$config->{JOB_FOLDER}/$job_id") or print_bug("Couldn't move to job folder for '$job_id': $!");

	$complete = -e "finished.flag";
    }

    unless ($complete) {
	# We can't check the pid until we fix how the forking is working
	# We'll settle for seeing if the error file is empty instead.

	if (-e "$config->{JOB_FOLDER}/$job_id/errors.txt") {
	    if ((stat "$config->{JOB_FOLDER}/$job_id/errors.txt")[7]) {
		print_bug("Job $job_id generated errors");
	    }
	}
	# We can check to see that the pid for this process is still alive
#	if (-e "$config->{JOB_FOLDER}/$job_id/pid.txt") {
#	    open(PID,"$config->{JOB_FOLDER}/$job_id/pid.txt") or print_bug("Couldn't open pid file for $job_id: $!");
#	    my $pid = <PID>;
#	    close PID;

#	    unless (kill 0, $pid) {
#		print_bug("Rscript for job $job_id died prematurely");
#	    }
#	}
    }

    my $template = HTML::Template -> new(filename => "$RealBin/../templates/results.html");

    $template -> param(JOB_ID => $job_id,
		       EXISTS => $exists,
		       COMPLETE => $complete);

    $template -> param(VERSION => $config->{VERSION},
		       ADMIN_EMAIL => $config->{ADMIN_EMAIL},
		       ADMIN_NAME => $config->{ADMIN_NAME});


    if ($complete) {
	# We can collect the data from the results folder

	# Read the hit table
	my @hit_table;
	open(IN,"GO_analysis_results.txt") or print_bug("Couldn't open GO results for $job_id: $!");
	$_ = <IN>; # Remove header

	while (<IN>) {
	    chomp;
	    my ($go_name,$query_count,$background_count,$category_count,$enrichment, $pval,$fdr,$potential_bias) = split(/\t/);

	    $go_name =~ s/\%/ /g;
	    
	    push @hit_table, {
		GO_NAME => $go_name,
		QUERY_COUNT => $query_count,
		BACKGROUND_COUNT => $background_count,
		CATEGORY_COUNT => $category_count,
		FDR => $fdr,
		POTENTIAL_BIAS => $potential_bias,
		ENRICHMENT => $enrichment,
	    };

	}

	$template -> param(HIT_TABLE => \@hit_table);
							


    }


    print $template -> output();

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

sub read_config {

    # Eventually we'll read this from our conf file, but let's hard code some stuff for now


    my $config = {

	ADMIN_EMAIL => 'simon.andrews@babraham.ac.uk',
	ADMIN_NAME => 'Simon Andrews',
	VERSION => '0.1.devel',
	JOB_FOLDER => "$RealBin/../../jobs/",
    }


}

