#!/usr/bin/perl
use warnings;
use strict;
use CGI;
use FindBin qw($RealBin);
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
    
    print $template -> output();


}


sub show_job {

}


