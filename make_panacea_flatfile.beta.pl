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

make_panacea_flatfile.pl 

=head1 SYNOPSIS

make_panacea_flatfile.pl [options] -d directory -t type -o outputFlatFile

	Options:
		-i	--input 	Input directory of the Pangenome information. Required.
		-o	--output 	Output file. Default is PanACEA.[date].txt
		-t	--type		Input type of the Pangenome. Currently supports "Pangenome" <default>, Pangenome "Iterative", "Panoct", and "Roary".
	
	
=head1 DESCRIPTION
make_panacea_flatfile.pl generates the panchromosome PanACEA flat file from pangenome output 
from  various programs, including PanOCT and Roary. This file is required to run the PanACAExT 
visualizer.

=head1 AUTHOR 
Thomas Clarke (tclarke@jcvi.org)
 2017 J. Craig Venter Institute

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

#load information from the files

my %types = ( "pangenome" => 1, "iterative" => 1, "panoct" => 1, "roary" => 1 );

my $core_cnt   = 0;
my $CORE_PERCT = 0.75;

my %gene_num;    #Matches the gene cluster name with the corresponding array number in the JSON matrix

my %clusters;    #stores names and number of genomes in cluster

my %fgi_grps;         #type = circular/linear; bound = the start and stop core genes; parent; sz = size; hit
my %fgi_member;       #the fGI that a cluster belongs to
my %cores;            #cores and their type
my %list;             #used to make fGI gene array: contains bounding genes, and the different fGRs
my %core_list;        #list of core genes
my %fgi_list;         #c_start, c_end
my %core_size;        #gets the core size for the number
my %fgi_gene_list;    #hash of arrays
my %geneLenInfo;      # hash containing the gene length mean, min, max and standard deviation for each cluster;
my %seqs;

#Directory to write the output
my $in_dir;           #Default value for the output head directory

my $date = time();

my $output_id = "PanACEA.txt";    #Default value for the output file
my $type      = "pangenome";
my $verb      = 1;
my $use_genbank	= 0;

#Doing the options here: need in_dir, out_dir, dir, func_file, output_id
my $help;
GetOptions(
    "input|i:s"  => \$in_dir,
    "type|t:s"   => \$type,
    "output|o:s" => \$output_id,
    "genbank|g"  => \$use_genbank,
	"help|h|?"   => \$help,
	
) or die "Can't get options!\n";

#Probably should set up a verbose tag to level the printing
if ( !$in_dir ) {

    warn "Please enter in a directory from which to create the Flat File\n";
    pod2usage( { -verbose => 2, -exitval => 0 } );

}

if ( !-e $in_dir ) {

    warn "$in_dir does not exist. Please check...\n";
    pod2usage( { -verbose => 2, -exitval => 0 } );

}

if ( !-d $in_dir ) {

    warn "$in_dir is not a directory. Please check...\n";
    pod2usage( { -verbose => 2, -exitval => 0 } );

}

if ( !$types{ lc($type) } ) {

    my $tmp = "$type is not a recognized pangenome type. Please check the below and select the correct one.\n";
    $tmp .= "Pangenome outputs able to be used: ";
    $tmp .= join( ", ", keys %types );
    die "$tmp\n";

}

if ( $verb == 1 ) {

    warn "Welcome to PanACEA Creation Software\n\n";
    warn "Using pangenome information located in $in_dir...\n\n\n";

}

if ( $output_id eq "PanACEA.txt" ) {

    if ( -e $output_id ) {

        $output_id = "PanACEA.$date.txt";
        warn "PanACEA.txt already exists. Using $output_id instead...\n\n";

    }

} else {

    if ( -e $output_id ) {

        warn "Warning: $output_id already exists, so it will be re-written...";

    }

}

if ( lc($type) eq "pangenome" ) { read_in_pangenome(); }
if ( lc($type) eq "iterative" ) { read_in_pangenome(); }
if ( lc($type) eq "panoct" )    { read_in_panoct(); }
if ( lc($type) eq "roary" )     { read_in_roary(); }

write_files();

if ( $verb == 1 ) {
    warn "Finished...\n\n";
}

exit(0);

sub read_in_panoct {

    warn "PanOCT function is not currently working. Please check again later..\n";

    #Pass data to gene_order
    exit(1);

}

#Reading in Roary output
sub read_in_roary {

    warn "Roary function is not currently working. Please check again later..\n";

    #load in x,y,z

    #make: 	cluster_wieghts
    #		num_core_adjacency_vector.txt consensus
    #		num_core_adjacency_
    exit(1);

}

#Reading in the Pangenome pipeline
sub read_in_pangenome {

    my ( $file_in1, $file_in2, $conses_in, $fgi_in, $fgi_arr, $comb_file, $frame_shift_file, $fasta_file, $sing_file );

    if ( -e $in_dir ) {

        if ( -d $in_dir ) {

            my $res_dir = "results_all";
            if ( !-e "$in_dir/$res_dir/" ) { $res_dir = "results"; }

            my $fgi_dir = "fGI";
            if ( !-e "$in_dir/$res_dir/$fgi_dir" ) { $fgi_dir = "fGIs"; }

            $file_in1  = "$in_dir/$res_dir/$fgi_dir/Core.attfGI";            #
            $file_in2  = "$in_dir/$res_dir/shared_clusters.txt";             #
            $conses_in = "$in_dir/$res_dir/$fgi_dir/consensus.txt";          #
            $fgi_in    = "$in_dir/$res_dir/$fgi_dir/FGI_inserts.details";    #
            if ( !-e $fgi_in ) { $fgi_in = "$in_dir/$res_dir/$fgi_dir/fGI_report.txt.details"; }    # ./results/fGI/

            $fgi_arr          = "$in_dir/$res_dir/$fgi_dir/fGI.att";          # ./results/fGI/ #
            $comb_file        = "$in_dir/combined.att";                       #./ Connects genes to different assemblies
            $frame_shift_file = "$in_dir/$res_dir/frameshifts.txt";           #./results/ For when genes get rename
            $sing_file        = "$in_dir/$res_dir/singletons_clusters.txt";

        } else {
            die "$in_dir is not a directory. Please double check...\n\n";
        }

    } else {
        die "Cannot find $in_dir for reading in pangenome files. Please double check...\n\n";
    }

    if ($verb) { warn "Reading in PanGenome information...\n"; }

    $in_dir =~ s/([\/]+)\Z//g;
    open( CORE1, "<", $file_in1 )
      or die("Cannot read in the Core + fGI attributes file. Please check $file_in1.\n");    #Core+fGI attributes file
    open( CLUSTER, "<", $file_in2 )
      or die("Cannot read in the Cluster file. Please check $file_in2.\n");                  #Shared Cluster information file
    open( SINGLETON, "<", $sing_file ) or die("Cannot read in the Singleton Gene file. Please check $sing_file.\n");
    open( CONSENS,   "<", $conses_in ) or die("Cannot read in the Consensus file. Please check $conses_in.\n");        #
    open( COMB,      "<", $comb_file ) or die("Cannot read in the Combined File file. Please check $comb_file.\n");
    open( FGI, "<", $fgi_in ) or die("Cannot read in the FGI details file. Please check $fgi_in.\n");    #Details file
    open( CORE2, "<", $fgi_arr )
      or die("Cannot read in the FGI Attributes file. Please check $fgi_arr.\n");                        #fGI attributes file
    open( FRAME_SHIFT, "<", $frame_shift_file )
      or die("Cannot read in the Frameshift file. Please check $frame_shift_file.\n");

    #get all the old sequence names in cases of frameshifts
    #oldGeneName is the hash we used to store the old name -> new name
    my %oldGeneName;

    while ( !eof(FRAME_SHIFT) ) {

        my $in = <FRAME_SHIFT>;
        chomp($in);
        my @in = split "\t", $in;
        for ( my $i = 1 ; $i < scalar(@in) ; $i++ ) {

            $oldGeneName{ $in[0] } = $in[$i];
            $oldGeneName{ $in[$i] } = $in[0];

        }

    }

    my %gene2clust;    #Hash with keys (1) gene name and (2) genomes, it gives value the cluster number
    my %genomes;

    warn "Reading in Clusters...\n";

#reads in the cluster id and the gene names associated with that cluster. From that gets the fasta sequenece associated with each cluster
#and with each cluster also the mean, sd, min and maxv length (it does take a while to run-> is that the calculation or in loading the multi-fasta?)
    my $in = <CLUSTER>;
    chomp($in);
    my @head = split "\t", $in;
    for ( my $i = 1 ; $i < scalar(@head) ; $i++ ) {

        $head[$i] =~ s/\s/\_/g;
        $genomes{ $head[$i] } = $i;

    }
    my $c = 0;
    while ( !eof(CLUSTER) ) {

        my $in = <CLUSTER>;
        chomp($in);
        my @n = split "\t", $in;
        for ( my $i = 1 ; $i < scalar(@head) ; $i++ ) {

            if ( $n[$i] ) {

                $n[$i] =~ s/'//g;
                $n[$i] =~ s/`//g;
                $n[$i] =~ s/\"//g;

                $clusters{ $n[0] }->{ $head[$i] } =
                  $n[$i];   #for each cluster id and genome returns the gene name; also includes centroid and number_of_genes

            }

        }
        $c++;
        $n[2] =~ s/\;//g;
        $n[2] =~ s/\'//g;

        #Fills  up the cluster ids per genome/gene name; also gets the
        my $minLen = -1;
        my $maxLen = -1;
        my $sum_len;
        my $sum_len2;
        my $n;

        for ( my $i = 7 ; $i < scalar(@head) ; $i++ ) {

            if ( $n[$i] ) {

                $gene2clust{ $head[$i] }->{ $n[$i] }       = $n[0];
                $seqs{byClust}->{ $n[0] }->{ $head[$i] }   = 1;
                $seqs{byClustId}->{ $n[0] }->{ $head[$i] } = $n[$i];

            }

        }
        $gene_num{ $n[0] } = $c;

        #adding
        #fills out the gene jso array: for this one gene cluster equals one row
        # columns are 1: ID; 2: number; 3: 4: 5: function

        #Tries to get the number of genomes in the total. Probably should be the scalar(@head)-7;
    }

    warn "Finished Reading in Clusters...\n";
    warn "$c genes clusters identified...\n";

    #Also get the ID of the singletons to connect geneID to cluster ID;
    $in = <SINGLETON>;
    chomp($in);
    my @head2 = split "\t", $in;
    for ( my $i = 1 ; $i < scalar(@head2) ; $i++ ) {

        $head2[$i] =~ s/\s/\_/g;

    }
    while ( !eof(SINGLETON) ) {

        my $in = <SINGLETON>;
        chomp($in);
        my @in = split "\t", $in;
        $gene2clust{ $in[6] }->{ $in[5] }        = $in[0];
        $clusters{ $in[0] }->{ $in[6] }          = $in[5];
        $seqs{byClust}->{ $in[0] }->{ $in[6] }   = 1;
        $seqs{byClustId}->{ $in[0] }->{ $in[6] } = $in[5];
        $clusters{ $in[0] }->{"centroid_locus"}  = $in[5];
        for ( my $i = 1 ; $i < scalar(@head2) ; $i++ ) {

            $in[$i] =~ s/'//g;
            $in[$i] =~ s/`//g;
            $in[$i] =~ s/\"//g;

            $clusters{ $in[0] }->{ $head2[$i] } =
              $in[$i];    #for each cluster id and genome returns the gene name; also includes centroid and number_of_genes

        }

    }

    my %genomeAttSort;    #hash with key as the contif
    my %genomeCnt;
    my %genedir;
	my %geneHits;
	warn "Reading and Sorting Genes...\n";

	if ( $use_genbank )
	{
		my $gb_dir = "$in_dir/gb_dir/";
		if (-d $gb_dir)
		{
			opendir(GBDIR, $gb_dir);
			my @gb_files = readdir(GBDIR);
			foreach my $gb_file (@gb_files)
			{
				if ($gb_file =~ /\A(\S+).gb\Z/)
				{
					my $genome_id 	= $1;
					my $st 			= "";
					my $end 		= "";
					my $gene_id 	= "";
					my $id        	= "";
					my $gene_prod 	= "";
					my $gene_dir  	= 1;
					warn "Opening file $gb_dir/$gb_file\n";
					open(GBFILE, "<", $gb_dir . "/" . $gb_file);
					my $var_cnt = 0;
					while (my $v = <GBFILE>)
					{
						$var_cnt++;
						if ($v =~ /\ALOCUS(\s+)(\S+)/)
						{
							my $locusID = $2; 
							if ($locusID =~ /\A(\S+)\_(\d+)\-(\d+)\Z/)
							{
								$st = $2;
								$end = $3;
								$id = $1;
							}
						}
						if ($v =~ /CDS(\t)complement/)
						{
							$gene_dir = -1;
						}
						if ($v =~ /locus\_tag\=\"([^\"]+)\"/)
						{
							$gene_id = $1;
						}
						if ($v =~ /\/product\=\"([^\"]+)\"/)
						{
							$gene_prod = $1;
						}
						if ($v =~ /\/\//)
						{
							if ($id && $st && $end)
							{
								$genomeCnt{$genome_id}++;
								if (ref($id) eq "HASH") {die("$id  $gene_id"); }
								$genomeAttSort{ $id }->{ min( $st, $end ) } = [ $id, $gene_id, $st, $end, $gene_prod, $genome_id, $gene_dir ];
								$geneHits{$genome_id}->{$gene_id} = $id;
							}
							$st 		= "";
							$end 		= "";
							$gene_id 	= "";
							$id        	= "";
							$gene_prod 	= "";
							$gene_dir  	= 1;
						}
					}
					close(GBFILE);
					warn "Read in " . $genomeCnt{$genome_id} . " loci in $genome_id\n";
				}
			}
		}
		else
		{
			warn "Cannot find genbank directory $gb_dir. Reverting to combined attributes file...\n";
			$use_genbank = 0;
		}
	}
	if ( !$use_genbank )
	{
		#reads throught the combined attributes files; doesn't have to be sorted as the script does that below
		while ( !eof(COMB) ) {

			my $v = <COMB>;
			chomp($v);
			my @n         = split "\t", $v;
			my $genome_id = $n[5];
			my $st        = $n[2];
			my $end       = $n[3];
			my $gene_id   = $n[1];
			my $id        = $n[0];
			my $gene_dir  = 1;
			if ( $st > $end ) { $n[2] = $end; $n[3] = $st; $gene_dir = -1; }
			$genomeCnt{$genome_id}++;
			$genomeAttSort{ $n[0] }->{ min( $st, $end ) } = [ $n[0], $n[1], $n[2], $n[3], $n[4], $n[5], $gene_dir ];

		}
	}
	
	foreach my $genomeIDs (keys(%gene2clust))
	{
		foreach my $geneIDs (keys(%{$gene2clust{$genomeIDs}}))
		{
			if (!$geneHits{$genomeIDs}->{$geneIDs})
			{
				#die("Cannot find $geneIDs from $genomeIDs\n");
			}
			else
			{
				warn("Did find $geneIDs from $genomeIDs\n");
			}
		}
	}
	
	
    my %clust2genome;    #Hash
    my %genome_order;    #hash of array
    my %fgi_2_genome;
	warn "Found ", scalar(keys(%genomeAttSort)), " contigs..\n";
	my @missingCount;
    foreach my $keys (keys(%genomeAttSort)) {

        #Sorts the
     	my @keys = sort { $b <=> $a } keys( %{ $genomeAttSort{$keys} } );
        for ( my $i = 0 ; $i < scalar(@keys) ; $i++ ) {

            my @n         = @{ $genomeAttSort{$keys}->{ $keys[$i] } };
            my $genome_id = $n[5];
            my $st        = $n[2];
            my $end       = $n[3];
            my $gene_id   = $n[1];
            my $id        = $n[0];
            my $gene_dir  = $n[6];
            if ( $genomes{ $n[5] } ) {

                if ( $gene2clust{$genome_id}->{$gene_id} ) {
					$missingCount[0]++;

                    $clust2genome{ "CL_" . $gene2clust{ $n[5] }->{$gene_id} }->{$genome_id} = "$st-$end";
                    if ( $genome_order{$genome_id}->{arr} ) {

                        push @{ $genome_order{$genome_id}->{arr} }, "CL_" . $gene2clust{$genome_id}->{ $n[1] };

                    } else {

                        $genome_order{$genome_id}->{arr}->[0] = "CL_" . $gene2clust{$genome_id}->{$gene_id};

                    }
                    $genome_order{$genome_id}->{list}->{ "CL_" . $gene2clust{$genome_id}->{$gene_id} } =
                      scalar( @{ $genome_order{$genome_id}->{arr} } );
                    $genome_order{$genome_id}->{id}->{ "CL_" . $gene2clust{$genome_id}->{$gene_id} } = $id;
                    $genedir{$genome_id}->{ $gene2clust{ $n[5] }->{$gene_id} } = $gene_dir;

                } else {

                    if ( $oldGeneName{$gene_id} ) {

                        $gene_id = $oldGeneName{$gene_id};
						
                    }
                    if ( $gene2clust{$genome_id}->{$gene_id} ) {
						$missingCount[0]++;
                        $clust2genome{ "CL_" . $gene2clust{ $n[5] }->{$gene_id} }->{$genome_id} = "$st-$end";
                        if ( $genome_order{$genome_id}->{arr} ) {

                            push @{ $genome_order{$genome_id}->{arr} }, "CL_" . $gene2clust{$genome_id}->{$gene_id};

                        } else {

                            $genome_order{$genome_id}->{arr}->[0] = "CL_" . $gene2clust{$genome_id}->{$gene_id};

                        }
                        $genome_order{$genome_id}->{list}->{ "CL_" . $gene2clust{$genome_id}->{$gene_id} } =
                          scalar( @{ $genome_order{$genome_id}->{arr} } );
                        $genome_order{$genome_id}->{id}->{ "CL_" . $gene2clust{$genome_id}->{$gene_id} } = $id;
                        $genedir{$genome_id}->{ $gene2clust{$genome_id}->{$gene_id} } = $gene_dir;

                    } else {
						$missingCount[1]++;
                        #print( STDERR "No hits..." . $n[5] . " " . $gene_id . " $id $st $end\n" );

                    }

                }

            } else {

                die("Cannot find genome " . $n[5] . " ... quitting\n\n");

            }

        }

    }
	warn "Found ", $missingCount[0], " genes and missing ", $missingCount[1], "...\n";
    warn "Reading in Core and FGI Cluster Assignments...\n";

    while ( !eof(CONSENS) ) {

        my $in = <CONSENS>;
        if ( $in =~ /\#Assembly_Core(\s+)(\d+)(\s+)(\S+)/ ) {

            $cores{$2}->{type} = $4;
            $core_cnt++;

        }

        if ( $in =~ /\#Assembly_fGI/ ) {

            chomp($in);
            my ( $type, $id_num, $dist, $gene_num, $nt_len ) = split "\t", $in;

            $fgi_grps{$id_num}->{type} = $dist;
            $in = <CONSENS>;
            my $st;
            my $end;
            my $hit;
            while ( $in =~ /\A[^\s\#]/ && !eof(CONSENS) ) {

                chomp($in);
                my @tmp = split "\t", $in;
                if ( $tmp[2] && $tmp[2] eq "P" ) {

                    $fgi_member{ $tmp[0] } = $id_num;
                    foreach my $genome ( keys( %{ $gene2clust{ $tmp[0] } } ) ) {

                        $fgi_2_genome{$id_num}->{$genome}++;

                    }
                    $hit = 1;

                }
                $in = <CONSENS>;

            }

        }

    }

    warn "Assigning genomes to FGR orders...\n";

    my %genome_fgi_list;
    my %genome_fgi_by_start
      ; #This list the FGI by starting Core Cluster. Contains keys: list (non-core genes in between); break (is there a break?)
        #st: the starting core cluster; end: the ending core cluster;
    foreach my $genome ( keys(%genome_order) ) {

        my $curr;
        my $tmp     = "";
        my $back    = 0;
        my $curr_id = "";
        my $tmp2    = "";
        my $back2   = 0;
        my $st_gene;
        my $st_gene_rev;
        my $st_gene_id;
        my $break;
        my $in_cnt;
		print STDERR"Running Genome: $genome\n";
        foreach my $gene_cl ( @{ $genome_order{$genome}->{arr} } ) {

            if ( $gene_cl =~ /CL_(\d+)/ ) {

                my $gene    = $1;
                my $geneDir = $genedir{$genome}->{$gene};
                my $num     = 3;
                my $revNum  = 5;
                if ( $geneDir == 1 ) { $num = 5; $revNum = 3; }

                #$revNum=$num;
                if ( $curr_id eq "" ) {

                    $curr_id = $genome_order{$genome}->{id}->{$gene_cl};
                    $break   = 1;

                }
                if ( $curr_id ne $genome_order{$genome}->{id}->{$gene_cl} && $st_gene ) {

                    if ($tmp) {

                        $genome_fgi_list{$tmp}->{$genome}                            = 1;
                        $genome_fgi_list{$back}->{$genome}                           = 1;
                        $genome_fgi_list{$tmp2}->{$genome}                           = 1;
                        $genome_fgi_list{$back2}->{$genome}                          = 1;
                        $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{list}  = $tmp2;
                        $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{break} = 1;
                        $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{st}    = $st_gene_id;

                        $genome_fgi_by_start{ "NA_" . $st_gene_rev }->{$genome}->{list}  = $back;
                        $genome_fgi_by_start{ "NA_" . $st_gene_rev }->{$genome}->{end}   = $st_gene_id;
                        $genome_fgi_by_start{ "NA_" . $st_gene_rev }->{$genome}->{break} = 1;

                    } else {

                        $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{list}  = "";
                        $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{break} = 1;

                        $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{st}       = $st_gene_id;
                        $genome_fgi_by_start{ "NA_" . $st_gene_rev }->{$genome}->{end}  = $st_gene_id;
                        $genome_fgi_by_start{ "NA_" . $st_gene_rev }->{$genome}->{list} = "";

                        $genome_fgi_by_start{ "NA_" . $st_gene_rev }->{$genome}->{break} = 1;

                    }
                    $back    = "";
                    $back2   = "";
                    $break   = 1;
                    $tmp     = "";
                    $tmp2    = "";
                    $st_gene = "";
                    $in_cnt  = 0;
                    if ( !$fgi_member{$gene_cl} ) {

                        $st_gene_id  = $gene;
                        $st_gene     = "$gene-$num";
                        $st_gene_rev = "$gene-$revNum";

                    }

                }
                $curr_id = $genome_order{$genome}->{id}->{$gene_cl};
                $curr    = $gene;
                #if ($st_gene_id ) { print STDERR ($st_gene . "_" . $gene . "-" . $revNum . " $genome\n"); }
				if ( !$fgi_member{$gene_cl} ) {

                    if ($tmp) {

                        $genome_fgi_list{ $tmp . $gene . ":" }->{$genome}  = 1;
                        $genome_fgi_list{ $tmp2 . $gene . ":" }->{$genome} = 1;
                        if ( $st_gene && ( $st_gene ne $gene . "-" . $revNum ) ) {

                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{list}  = $tmp2;
                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{st}    = $st_gene_id;
                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{end}   = $gene;
                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{break} = $break;

                            $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{list} = $tmp2;
                            $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{st}   = $st_gene_id;
                            $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{end}  = $gene;

                            $genome_fgi_by_start{ $st_gene . "_NA" }->{$genome}->{break} = $break;

                            $genome_fgi_by_start{ "NA_" . $gene . "-" . $revNum }->{$genome}->{list}  = $tmp2;
                            $genome_fgi_by_start{ "NA_" . $gene . "-" . $revNum }->{$genome}->{st}    = $st_gene_id;
                            $genome_fgi_by_start{ "NA_" . $gene . "-" . $revNum }->{$genome}->{end}   = $gene;
                            $genome_fgi_by_start{ "NA_" . $gene . "-" . $revNum }->{$genome}->{break} = $break;

                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{list} = $back;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{st}   = $gene;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{end}  = $st_gene_id;

                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{break} = $break;

                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_NA" }->{$genome}->{list} = $back;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_NA" }->{$genome}->{st}   = $gene;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_NA" }->{$genome}->{end}  = $st_gene_id;

                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_NA" }->{$genome}->{break} = $break;

                            $genome_fgi_by_start{ "NA_" . $st_gene }->{$genome}->{list} = $back;
                            $genome_fgi_by_start{ "NA_" . $st_gene }->{$genome}->{st}   = $gene;
                            $genome_fgi_by_start{ "NA_" . $st_gene }->{$genome}->{end}  = $st_gene_id;

                            $genome_fgi_by_start{ "NA_" . $st_gene }->{$genome}->{break} = $break;

                        } else {

                            $genome_fgi_by_start{ "NA_" . $gene . "-" . $revNum }->{$genome}->{list}  = $tmp2;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_NA" }->{$genome}->{list}  = $back;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_NA" }->{$genome}->{st}    = $gene;
                            $genome_fgi_by_start{ "NA_" . $gene . "-" . $revNum }->{$genome}->{end}   = $gene;
                            $genome_fgi_by_start{ "NA_" . $gene . "-" . $revNum }->{$genome}->{break} = $break;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_NA" }->{$genome}->{break} = $break;

                        }
                        $genome_fgi_list{$back2}->{$genome} = 1;
                        $genome_fgi_list{$back}->{$genome}  = 1;
                        $back                               = "";
                        $back2                              = "";
                        $tmp2                               = "";
                        $tmp                                = "";
                        $break                              = 0;
                        $in_cnt                             = 0;
                        $st_gene_id                         = $gene;
                        $st_gene_rev                        = $gene . "-" . $revNum;
                        $st_gene                            = $gene . "-" . $num;

                    } else {

                        if ( $st_gene && ( $st_gene ne $gene . "-" . $revNum ) ) {

                            $genome_fgi_list{ $gene . ":" }->{$genome}                                         = 1;
                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{list}  = "";
                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{end}   = $gene;
                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{st}    = $st_gene_id;
                            $genome_fgi_by_start{ $st_gene . "_" . $gene . "-" . $revNum }->{$genome}->{break} = $break;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{list}  = "";
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{end}   = $gene;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{st}    = $st_gene_id;
                            $genome_fgi_by_start{ $gene . "-" . $revNum . "_" . $st_gene }->{$genome}->{break} = $break;

                            $back        = "";
                            $back2       = "";
                            $tmp2        = "";
                            $tmp         = "";
                            $break       = 0;
                            $in_cnt      = 0;
                            $st_gene_id  = $gene;
                            $st_gene_rev = $gene . "-" . $revNum;
                            $st_gene     = $gene . "-" . $num;

                        }

                    }

                } else {

                    my $fgi_tmp = $fgi_member{$gene_cl};
                    $in_cnt++;
                    if ( !$tmp2 || $tmp2 !~ /$gene:/ ) {

                        $tmp2 .= "$gene:";
                        $back = "$gene:$back";

                    }
                    if ( $tmp !~ /$fgi_tmp/ ) {
                        $tmp .= "$fgi_tmp:";
                    }

                }

            } else {
                die("Cannot find any gene clusters under the ID $gene_cl... Quiting...\n");
            }

        }

    }

    warn "Reading in Core and FGI Cluster Assignments II...\n";

    my ( $prev_n, $prev );
    while ( !eof(CORE1) ) {

        my $in = <CORE1>;
        chomp($in);
        my @n = split "\t", $in;

        #$fgi_member{$n[1]} = $n[0];
        push( @{ $core_list{ $n[0] } }, \@n );
        $core_size{ $n[1] } = $n[3] - $n[2];

        if ( !$cores{ $n[0] }->{sz} ) {
            $cores{ $n[0] }->{sz} = $n[3];
        }

        if ( $n[3] > $cores{ $n[0] }->{sz} ) {
            $cores{ $n[0] }->{sz} = $n[3];
        }

        if ( $n[2] > $cores{ $n[0] }->{sz} ) {
            $cores{ $n[0] }->{sz} = $n[2];
        }

        if ( $prev && $prev =~ /INS/ ) {
            $fgi_list{$prev}->{c_end} = $n[1];
        }

        if ( $n[4] && $n[4] =~ /fGI/ ) {
            $fgi_list{ $n[1] }->{c_start} = $prev;
        }

        $prev = $n[1];

    }

    my ( $st, $prev2, $end, $hit );
    while ( !eof(CORE2) ) {

        my $in = <CORE2>;
        chomp($in);
        my @n = split "\t", $in;
        $fgi_grps{ $n[0] }->{sz} = 0;
        if ( $prev2 && $prev2 ne $n[0] ) {

            $fgi_grps{$prev2}->{bound} = $st . "_" . $end;
            my $st_1 = "";
            if ( $fgi_member{ "CL_" . $st } ) { $st_1 = $fgi_member{ "CL_" . $st }; }
            my $end_1 = "";
            if ( $fgi_member{ "CL_" . $end } ) { $end_1 = $fgi_member{ "CL_" . $end }; }

            $fgi_grps{$prev2}->{parent} = $st_1 . "_" . $end_1;
            $st                         = "";
            $end                        = "";
            $hit                        = "";

        }
        $prev2 = $n[0];
        if ( $n[1] !~ /CONTEXT/ ) {

            push( @{ $core_list{ $n[0] } }, \@n );
            $core_size{ $n[1] }     = $n[3] - $n[2];
            $fgi_gene_list{ $n[1] } = \@n;
            if ( $n[3] > $fgi_grps{ $n[0] }->{sz} ) {

                $fgi_grps{ $n[0] }->{sz} = $n[3];

            }
            if ( $n[2] > $fgi_grps{ $n[0] }->{sz} ) {

                $fgi_grps{ $n[0] }->{sz} = $n[2];

            }
            $prev = $n[1];
            $hit  = 1;

        } else {

            if ( !$hit && $n[1] =~ /CL_(\d+)/ ) { $st = $1; }
            if ( $hit && !$end && $n[1] =~ /CL_(\d+)/ ) { $end = $1; }

        }

    }

    my $ins_id;
    my $st_1   = "";
    my $end_1  = "";
    my $st_id  = "";
    my $end_id = "";
    open( FGI, "<", $fgi_in ) or die();
    while ( !eof(FGI) ) {

        my $in = <FGI>;
        chomp($in);
        my @n = split "\t", $in;
        my $line = $in;
        if ( $n[0] =~ /INS/ ) {

            my $cnt_i;
            $ins_id = $n[0];
            $st_id  = "";
            $end_id = "";

            if ( scalar(@n) > 1 ) {

                $n[1] =~ /(\d+)\_(\d+)/;
                $st_1  = $1 . "-" . $2;
                $st_id = $1;

            }
            if ( scalar(@n) > 2 ) {

                $n[2] =~ /(\d+)\_(\d+)/;
                $end_1  = $1 . "-" . $2;
                $end_id = $1;

            }


            if ( $st_id eq $end_id ) { die( $in . " $st_id $end_id" ); }
            my %out_list;
            my %hits;
            my $hit_type = 0;
			if ($ins_id eq "CL_INS_673") { print( "Here: $st_1 $end_1\n"); }
            foreach my $a ( keys( %{ $genome_fgi_by_start{ $st_1 . "_" . $end_1 } } ) ) {
                my $in = $st_id . ":" . $genome_fgi_by_start{ $st_1 . "_" . $end_1 }->{$a}->{list} . $end_id . ":";
			if ($ins_id eq "CL_INS_673") { die( "Here: $st_1 $end_1 "); }

                if (!$end_id) {
                    $end_1  =~ /(\d+)\-(\d+)/; 
                    $end_id = $1;
                }
                if (!$st_id) {
                    $st_1  =~ /(\d+)\-(\d+)/;
                    $st_id = $1;
                }
                if ($st_id ne $end_id) {

                    #warn "$st_1\t$end_1\t$in\n";
                    $out_list{$in}->{st}   = "CL_" . $st_id;
                    $out_list{$in}->{list} = $genome_fgi_by_start{ $st_1 . "_" . $end_1 }->{$a}->{list};
                    $out_list{$in}->{en}   = "CL_" . $end_id;
                    $out_list{$in}->{cnt}++;
                    $out_list{$in}->{gen_list} .= $a . ":";
                    $hits{$a} = 1;
                    $hit_type = 1;
                }

            }
            foreach my $a ( keys( %{ $genome_fgi_by_start{ $st_1 . "_NA" } } ) ) {
							if ($ins_id eq "CL_INS_673") { print( "Here: $st_1 $end_1\n"); }

                if ( !$hits{$a} && $st_id ) {

                    my $in    = $st_id . ":" . $genome_fgi_by_start{ $st_1 . "_NA" }->{$a}->{list};
                    my $oth   = "";
                    my $oth_t = "";
                    if ( $genome_fgi_by_start{ $st_1 . "_NA" }->{$a}->{end}
                        && !$fgi_member{ "CL_" . $genome_fgi_by_start{ $st_1 . "_NA" }->{$a}->{end} } ) {

                        #$in .= $genome_fgi_by_start{$st_1 . "_NA"}->{$a}->{end} .":";
                        $out_list{$in}->{en} = "CL_" . $genome_fgi_by_start{ $st_1 . "_NA" }->{$a}->{end};

                    } else {

                        if (   $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}
                            && $genome_fgi_by_start{ $st_1 . "_" }->{$a}->{break} ) {

                            $in .= $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{list} . $end_id . ":";
                            $hits{$a} = 1;
                            $out_list{$in}->{en} = "CL_" . $end_id;
                            $oth   .= $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{list};
                            $oth_t .= $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{type} . ":";

                        }

                    }
                    $out_list{$in}->{st} = "CL_" . $st_id;

                    #$out_list{$in}->{en} = "";
                    $out_list{$in}->{cnt}++;
                    $out_list{$in}->{gen_list} .= $a . ":";
                    $out_list{$in}->{list} = $genome_fgi_by_start{ $st_1 . "_NA" }->{$a}->{list} . $oth;
                    $hit_type = 2;

                }

            }
            foreach my $a ( keys( %{ $genome_fgi_by_start{ "NA_" . $end_1 } } ) ) {
						if ($ins_id eq "CL_INS_673") { print( "Here2: $st_1 $end_1\n "); }


                if ( !$hits{$a} && $end_id ) {

                    my $in = $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{list} . "$end_id:";

                    if ( $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{st}
                        && !$fgi_member{ "CL_" . $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{st} } ) {

                        #$in = $genome_fgi_by_start{ "NA_" . $end_1}->{$a}->{st} .":$in";
                        $out_list{$in}->{st} = "CL_" . $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{st};

                    }

                    $out_list{$in}->{cnt}++;
                    $out_list{$in}->{gen_list} .= $a . ":";
                    $out_list{$in}->{en} = "CL_" . $end_id;

                    #$out_list{$in}->{st} = "";
                    $out_list{$in}->{list} = $genome_fgi_by_start{ "NA_" . $end_1 }->{$a}->{list};

                    $hit_type = 3;

                }

            }
			
            foreach my $a ( keys(%out_list) ) {

                if ( $a eq ":" ) {

                    die( $hit_type . " " . $ins_id . " S-" . $st_id . " E-" . $end_id . " " . $out_list{$a}->{gen_list} );

                }
                my @in = split ":", $out_list{$a}->{list};
                $cnt_i->{$ins_id} += $out_list{$a}->{cnt};
                $list{$ins_id}->{$a}->{arr} = \@in;
                $list{$ins_id}->{$a}->{cnt} = $out_list{$a}->{cnt};
                $list{$ins_id}->{$a}->{st}  = $out_list{$a}->{st};
                $list{$ins_id}->{$a}->{end} = $out_list{$a}->{en};
                $list{$ins_id}->{$a}->{gen} .= $out_list{$a}->{gen_list};

            }

        } else {

            my $start;
            my $st_dir;
            my $end;
            my $end_dir;
            my $fam = "";
            my @in;
            $n[0] .= ":";
            my $chain  = "";
            my $chain2 = "";

            #my $hits; my $curr;
            while ( $n[0] =~ /(\D*)(\d+)([\-\+])(\:*)/g ) {

                my $core_id = $1;
                my $id      = $2;
                my $dir_pl  = $3;
                my $dir;
                if   ( $dir_pl eq "+" ) { $dir = 1; }
                else                    { $dir = -1; }
                if ($core_id) {

                    $chain  .= "$id:";
                    $chain2 .= "$id:";
                    if ( $core_id =~ /START/ ) {

                        $start  = "CL_$id";
                        $st_dir = $dir;
                        if ( $id eq $end_1 ) {

                            $end     = "CL_$id";
                            $end_dir = 1;
                            if ( $dir_pl eq "+" ) { $end_dir = -1; $dir = 1; }
                            else                  { $dir = -1; }
                            $start  = "";
                            $st_dir = 0;

                        }

                    }

                    if ( $core_id =~ /STOP/ ) {

                        $end     = "CL_$id";
                        $end_dir = $dir;

                    }

                    if ( $core_id =~ /U/ ) {

                        $end     = "CL_$id";
                        $end_dir = $dir;

                    }

                } else {

                    if ( !$fgi_member{ "CL_" . $id } ) { die($id); }
                    my $tst = $fgi_member{ "CL_" . $id } . ":";
                    if ( $fam !~ /$tst/ ) { $fam .= $tst; $chain .= $tst; }
                    push @in, $id;
                    $chain2 .= "$id:";

                }

            }

        }

    }

}

sub write_files {

    open( OUT_FILE, ">", $output_id ) or die "Cannot open $output_id. Please check...\n";
    if ( $verb == 1 ) {
        warn "Writing PanACEA flat file to $output_id...\n";
    }

    #Writing each of the modules

    #Starting the chromosomes
    for my $cur_chr ( keys(%cores) ) {

        print OUT_FILE "START\tCHROM\t$cur_chr\n";
        print OUT_FILE "type\t", $cores{$cur_chr}->{type}, "\n";
        print OUT_FILE "sz\t",   $cores{$cur_chr}->{sz},   "\n";
        print OUT_FILE "is_core\t1\n";
        print OUT_FILE "END\n";

        #In each chromosome, making the core regions
        foreach my $core ( @{ $core_list{$cur_chr} } ) {

            my ( $chr, $ID, $st, $end, $def, $type, $len ) = @$core;
            print OUT_FILE "START\tCORE\t$ID\n";
            print OUT_FILE "chr\t$chr\n";
            print OUT_FILE "start\t$st\n";
            print OUT_FILE "end\t$end\n";
            print OUT_FILE "def\t$def\n";
            print OUT_FILE "type\t$type\n";
            print OUT_FILE "len\t$len\n";
            print OUT_FILE "size\t", $core_size{$ID}, "\n";
            print OUT_FILE "END\n";

        }

    }

    #Writing each of the fGR regions that are shown on the pan-chromosomes
    foreach my $id ( keys(%list) ) {

        foreach my $id2 ( keys( %{ $list{$id} } ) ) {

            print OUT_FILE "START\tFGR\t$id\n";
            print OUT_FILE "order\t", $id2, "\n";
            if ( $list{$id}->{$id2}->{st} ) {

                print OUT_FILE "st\t", $list{$id}->{$id2}->{st}, "\n";

            }
            if ( $list{$id}->{$id2}->{end} ) {

                print OUT_FILE "end\t", $list{$id}->{$id2}->{end}, "\n";

            }
            print OUT_FILE "gen\t", $list{$id}->{$id2}->{gen}, "\n";

            print OUT_FILE "END\n";

        }

    }

    #Writing the gene cluster modules
    foreach my $id ( keys(%clusters) ) {

        print OUT_FILE "START\tCLUSTER\t$id\n";

        print OUT_FILE "protein_name\t",   $clusters{$id}->{"protein_name"},   "\n";
        print OUT_FILE "num_of_members\t", $clusters{$id}->{"num_of_members"}, "\n";
        print OUT_FILE "centroid\t",       $clusters{$id}->{"centroid_locus"}, "\n";
        my $genome_list;
        my $name_list;
        my $dir_list;
        foreach my $gen_id ( keys( %{ $seqs{byClust}->{$id} } ) ) {

            $genome_list .= $gen_id . ";";
            $name_list   .= $seqs{byClustId}->{$id}->{$gen_id} . ";";

        }
        if ( $fgi_gene_list{ "CL_" . $id } ) {

            print OUT_FILE "start\t", $fgi_gene_list{ "CL_" . $id }->[2], "\n";
            print OUT_FILE "end\t",   $fgi_gene_list{ "CL_" . $id }->[3], "\n";

        }
        print OUT_FILE "genomes\t", $genome_list, "\n";
        print OUT_FILE "names\t",   $name_list,   "\n";
        print OUT_FILE "END\n";

    }

    #Writing the fGI clusters
    foreach my $id ( keys(%fgi_grps) ) {

        if ( $fgi_grps{$id}->{type} eq "cycle" ) {

            print OUT_FILE "START\tCHROM\t$id\n";
            print OUT_FILE "type\t", $fgi_grps{$id}->{type}, "\n";
            print OUT_FILE "sz\t",   $fgi_grps{$id}->{sz},   "\n";
            print OUT_FILE "is_core\t0\n";
            print OUT_FILE "END\n";
            foreach my $core ( @{ $core_list{$id} } ) {

                my ( $chr, $ID, $st, $end, $def, $type, $len ) = @$core;
                print OUT_FILE "START\tCORE\t$ID\n";
                print OUT_FILE "chr\t$chr\n";
                print OUT_FILE "start\t$st\n";
                print OUT_FILE "end\t$end\n";
                print OUT_FILE "def\t$def\n";
                print OUT_FILE "type\t$type\n";
                print OUT_FILE "len\t$len\n";
                print OUT_FILE "size\t", $core_size{$ID}, "\n";
                print OUT_FILE "END\n";

            }

        }

    }

}

sub min {

    if ( $_[0] > $_[1] ) { return $_[1]; }
    return ( $_[0] );

}

sub max {

    if ( $_[0] < $_[1] ) { return ( $_[1] ); }
    return ( $_[0] );

}

