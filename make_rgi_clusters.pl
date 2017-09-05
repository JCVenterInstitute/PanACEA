#!/usr/bin/env perl
#Copyright (C) 2016-2021 The J. Craig Venter Institute (JCVI).  All rights reserved

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>

=pod

=head1 NAME 

make_rgi_clusters.pl 

=head1 SYNOPSIS

make_rgi_clusters.pl [options] -d directory -t type -o outputFlatFile

	Options:
		-d	--input 	Location of dataSummary.txt file that is read. Required.
		-o	--output 	Output file. Default is aro_centroid.list.txt.
		-t	--type		Output type of the file. "Best" (all) or "all".
	
	
=head1 DESCRIPTION

make_rgi_clusters.pl generates the PanACEA readable term map file from a RGI generated
dataSummary.txt. Can make either a map film showing only the best hit (the default) or
all hits.

=head1 AUTHOR 

 Thomas Clarke (tclarke@jcvi.org)
 2017 J. Craig Venter Institute

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Usage perl make_rgi_clusters.pl dataSummary.txt output_centroid_list.txt
# output is
# centroid_id	AMR	ARO1,ARO2,etc

my $rgi_file;
my $output_file = "aro_centroid.list.txt";
my $type        = "best";
my $help;

GetOptions(
    "input|i:s"  => \$rgi_file,
    "type|t:s"   => \$type,
    "output|o:s" => \$output_file,
    "help|h|?"   => \$help,
) || die "Error getting options!\n";

pod2usage( { -verbose => 2, -exitval => 0 } ) if ( $help );

#Probably should set up a verbose tag to level the printing
if ( !$rgi_file ) {

    warn "Please enter in a directory from which to create the ARO map file\n";
    pod2usage( { -verbose => 2, -exitval => 0 } );

}

if ( !-e $rgi_file ) {

    warn "$rgi_file does not exist. Please check...\n";
    pod2usage( { -verbose => 2, -exitval => 0 } )

}

$type = lc($type);

if ( $type ne "best" && $type ne "all" ) {

    warn "Cannot recognize output type. Reverting to default (best)...\n\n";
    $type = "best";

}

my $list;

open( my $ifh, "<", $rgi_file )    or die("Cannot open rgi_file: $!\n");
open( my $ofh, ">", $output_file ) or die("Cannot open output file: $!\n");

while ( <$ifh> ) {

    next if /ORF_ID/;

    chomp;
    my @arr = split "\t", $_;
    if ( $arr[0] =~ /centroid_(\d+)/ ) {

        my $id    = "centroid_$1";
        my $best  = $arr[8];
        my @names = split ", ", $arr[11];
        my @aro   = split ", ", $arr[10];

        if ( $type eq "best" && $best ) {

            for ( my $i = 0 ; $i < scalar(@aro) ; $i++ ) {

                if ( $names[$i] && $names[$i] eq $best ) {

                    $aro[$i] =~ s/(\s)//g;
                    $list->{$id}->{ $aro[$i] }++;

                }

            }

        } else {

            while ( $arr[10] =~ /([^,]+)/g ) { 

                my $aro = $1; 
                $aro =~ s/(\s)//g; 
                $list->{$id}->{$aro}++; 

            }

        }

    }

}

foreach my $id ( keys(%$list) ) {

    my $out = join( ',', keys( %{ $list->{ $id } } ) );
    print $ofh "$id\t$out\n";

}
