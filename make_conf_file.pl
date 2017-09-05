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
use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path';
use Cwd 'cwd';

my ( $in_dir, $out_dir, $aro_loc, $go_loc, $help );
GetOptions( 'dir|d=s'  => \$in_dir, 
            'out|o=s'  => \$out_dir, 
            'go|g=s'   => \$go_loc,
            'aro|a=s'  => \$aro_loc,
            'help|h|?' => \$help );

if ( $help || !($in_dir) || !( -d $in_dir ) ) {

    my $HELPINFO=<<"EOF";
Runs MASH of a directory or fasta file against the MASH sketch\n
--------------------USAGE--------------------------------------
    -d input pangenome directory
    -o output directory
    -g location of go term obo file. Default is input directory
    -a location of aro obo file. Default is input directory
    -h help
EOF

    warn "$HELPINFO\n";
    exit(0);

}

my @cwd = split "\/", abs_path($0);

my $src_dir = "";
for ( my $i = 0 ; $i < scalar(@cwd) - 1 ; $i++ ) { $src_dir .= $cwd[$i] . "/"; }
chop $src_dir;

$in_dir = abs_path( $in_dir );

my $num       = 2;
my @obo_files = ( "$in_dir/aro.obo",               "$in_dir/obo/gene_ontology.1_2.obo" );
my @types     = ( "AROTerms",                      "GOTerms" );
my @map_files = ( "$in_dir/aro_centroid.list.txt", "$in_dir/centroids.mapterm.txt" );
my @color     = ( "ff0000",                        "$src_dir/cluster_colors.txt" );
my @names     = ( "Antibiotic Resistance",         "$in_dir/results/centroids.cluster_roles.txt" );

if ( !-e $names[1] ) {
    $names[1] = "DIR/centroids.cluster_roles.txt";
}

if ($aro_loc) {

    if ( -e $aro_loc ) {
        $obo_files[0] = $aro_loc;
    } else {
        warn "Cannot find aro obo file. Defaulting to input directory...";
    }

}

if ($go_loc) {

    if ( -e $go_loc ) {
        $obo_files[1] = $go_loc;
    } else {
        warn "Cannot find aro obo file. Defaulting to input directory...";
    }

}


if ( !$out_dir ) {
    $out_dir = cwd();
}

if ( !-e $out_dir ) {

    die "Cannot find directory. Quiting...\n\n\n";

}

if ( -e "$out_dir/func_file.conf.txt" ) {

    warn "Functional file already exists... overwriting...\n";

}

open( my $ofh, ">", "$out_dir/func_file.conf.txt" );

for ( my $i = 0 ; $i < $num ; $i++ ) {

    if ( $obo_files[$i] && $types[$i] && $map_files[$i] && $color[$i] && $names[$i] ) {

        printf $ofh "START\n";
        printf( $ofh "ID\t%s\n", $types[$i] );
        printf( $ofh "mapfile\t%s\n", $map_files[$i] );
        printf( $ofh "ontology\t%s\n", $obo_files[$i] );
        printf( $ofh "name\t%s\n", $names[$i] );
        printf( $ofh "color\t%s\n", $color[$i] );
        printf( $ofh "END\n\n" );

    }
    
}
