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

make_panacea.pl

=head1 SYNOPSIS

make_panacea.pl [options] -i PanACEA.InputFile.txt

	Options:
		-i	--input 	Input PanACEA Flat File. Required.
		-o	--output 	Output Directory for all the files. Default is current directory.
		-n	--name		Header name of HTML and SVG files. Default is \"PanHTML\"
		-d	--dir_name	Root directory of the HTML for use in the browser. Default is
					current directory.
		-f	--function	File containing information about additional gene cluster annotation.
		-a	--fasta		Either directory containing the cluster multiple alignment files
					or a single fasta file. Multiple inputs could be added in a comma separated
					list
		-t	--tree		File containing the phylogeny of the genomes in newick format.
					(Requires BioPerl to use)
		-m	--metafiles	File containing group metadata for the genomes. Only used in
					connection with tree.
		-g	--graphic	File containing the graphics configure file to change the
					size of the whole image and parts of the image


=head1 DESCRIPTION

make_panacea.pl generates an interconnected set of html, svg, javascript and json files
from a PanACEA flat file generated from a possible suite of programs, including Pangenome and
Roary, that allows users to explore multiple related prokayrotic genomes through their
pan-chromosomes, including core and flexible regions. While the basic version only shows
the location and relative position of the core and flexible regions, the user can also
add information on the gene cluster annotation and alignment as well as the genomic
annotation and phylogeny to obtain a more complete view of the relationship between the
genomes, their function and any regions with their associated genes that might drive
this relationship. The end result is a central html file that can be opened in many web
browsers.

=head2 Notes

=head3 'function' annotation config file format (-f flag)

PanACEA can use many types of annotations , but requires each to be inputed in modular format
in a single configure file with the format shown below. Each annotation block begins with
START and ends with END, and includes:
	* START ID	a string used on the table tab
	* mapfile:	a tab-delineated file with gene or centroid ID in the first column
			and a comma-separated term list in the second
	* ontology:	an obo type onotology file containing definitions for the terms in
			the mapfile
	* name:		a tab-delineated file with gene or centroid ID in the first column
			and a function in the second, OR a string that is assigned to each of the genes in the mapfile
	* color:	a tab-delineated file with the function in the first column and a color
			in the second, OR a rgb color string that is assigned to each of the gene
			with the functionq

More information can be found in the manual

=head3 'graphic' configure file (-g flag)

The graphical annotation config file allows the user to  change the dimensions and colors
of the image. The file consists of unique lines with the variable name followed by a space
and then followed by the value for one of 6 different variables. Variable names not listed
below will be ignored. Repeated variables will be given the last value. A list of the
variables along with a description follows:
	* BACKGROUND	color of the background of the outer ring/core regions on the main
					view
	* fGR			color of the fGRs in the outer ring/core regions on the main view
	* BORDER_SIZE	size of the border between the outer ring and the edge of the
					circle in pixels
	* CIRCLE_SIZE	size of the radius of the outer edge of the outer ring in pixels
	* TEXT_SIZE		size of the text in pixels
	* GENE_SIZE		size of the gene in pixels
	* DPI			DPI of the output PNGs

=head1 AUTHOR

 Thomas Clarke (tclarke@jcvi.org)
 2017 J. Craig Venter Institute

=cut

#Detailed decriptions of the sub-ruetines are found at the end of the script

use strict;
use Getopt::Long;
use Cwd;
use warnings;
use Pod::Usage;

#load information from the files

my $GENOME_NUM = 0;

#array function
my @funcInfo;

#This an image of the disk
#Used on a couple of images
my $disk_image_orig = disk_image("default");

#This is the image of the disk appropiate for javascript
my $disk_image_js = disk_image("js");

my $core_cnt   = 0;       #Running count of the number of core regions to make the figures
my $CORE_PERCT = 0.75;    #The core percentage- used to get a minimal y-axis in the core region graphs

my %jso;                  #hash used to make the tables to be made into JSON format for javascript to read
my %head;                 #keeps track of the table headers for each of the tables
my %col_width;            # the widths
my %gene_num;             #Matches the gene cluster name with the corresponding array number in the JSON matrix
my %seqs;                 #keeps the sequences associated with genes
my %gene_func;            #gene function
my %clusters;             #stores names and number of genomes in cluster
my %fgi_grps;             #type = circular/linear; bound = the start and stop core genes; parent; sz = size; hit
my %fgi_member;           #the fGI that a cluster belongs to
my %cores;                #cores and their type
my %list;                 #used to make fGI gene array: contains bounding genes, and the different fGRs
my $seq_len = 0;          #The length of genome
my %core_list;            #list of core genes
my %tick;                 #gets information for ticks around core
my %core_size;            #gets the core size for the number
my %fgi_gene_list;        #hash of arrays
my %geneLenInfo;          # hash containing the gene length mean, min, max and standard deviation for each cluster;
my %jsHead;               #json version of the table headers
my %go_name;              #names for the go terms (and any other added funnction) and the genes associated with it

my $tree_level = 0;       #Level at which trees are initial shown with 1 = leaves. O = tree not loaded
my %meta_types;           #Differnt genome metadata variables
my $tree_size = 500;      #Size of the tree image square
my %func_list;            #list of the functions associated with the genes
my $defaultFunct;         #The default funtional value

#### User Inputed variables on Command line###
my $mult_align_dir = "";           #The directory containing the cluster multiple alignments
my $dir            = cwd();        #Directory to write the output
my $out_dir        = cwd();        #Default value for the output head directory
my $output_id      = "PanACEA";    #Default value for -n
my $func_file;                     #Location of the configure file with the
my $file_in1;                      #Location of the PanACEA flat file
my $verb = 1;                      #Verbosity level: not currently used
my $graphic_config_file;           #Configure file with graphical output information
my $tree_file;                     #File with newick tree
my $meta_data_file;                #file(s) with group information with the newick trees

#Graphical values
# default values:
my $PI = 3.14159265358979323846;

my $max_radius   = 400;                 #Radius of the chromosome image
my $border       = 120;                 #Size of the border between chromosome and the edge of the screen
my $layer_height = $max_radius / 12;    #
my $gene_height  = 24;                  #height of the gene icons in the image
my $tick_height  = $gene_height / 3;    # tick size on the edge of the image
my $text_height  = 6;                   #Height of the
my $thumb_dif    = 0.15;                #The proportional size of the thumbnail sketches
my $DPI          = 3;                   #gives the DPI for outputed PNG images

my $SVGHEIGHT = 2 * ( $border + $max_radius );    #Height of the total image
my $SVGWIDTH  = 2 * ( $border + $max_radius );    #Width of the total image
my $thumb_size_h = $SVGHEIGHT * $thumb_dif;                 # size of the thumbnails
my $thumb_size_w = $SVGWIDTH * $thumb_dif;                  # size of thumbnails
my $chng_w       = $gene_height + $thumb_size_w;            #
my $disp         = "inline-block";                          #The
my $bw           = ( $max_radius * 2 - 2.25 * $border );    #body width
my $bodyTop      = "30px";                                  #Sets the location of the table body
my $bodyHeight   = ( $max_radius / 4 - 30 ) . "px";         #height
my $bodyWidth    = "100%";                                  #width of the table body

my $cr_start = 0;                                           #starting number of chromosomes
my $cr_cnt;                                                 #number of chromosomes

#Default name for files generated by the program
my $out_name    = "OUT";                                    #
my $filenameSVG = $out_name . ".svg";                       #written svgs
my $filenameTXT = $out_name . ".txt";                       #written txt files
my $filenamePNG = $out_name . ".png";                       #written png files

my %term2type;
my %type2term;
my $cur_st;  #current start value when drawing the pan chromosomes: this allows the horseshoe linear chromosomes to be drawn
my $default_table      = "ButtonIDRegion";    #The starting table selection Javascript ID
my $default_table_type = "Region";            #The starting table selection
my %ont;                                      #Connect term IDs to their OBO definitions
my %rgb;                                      #Connects terms IDS to their color

my $help;

#Getting all the user inputed usage
GetOptions(
    "i|input:s"    => \$file_in1,
    "t|tree:s"     => \$tree_file,
    "m|metadata:s" => \$meta_data_file,
    "a|fasta:s"    => \$mult_align_dir,
    "o|output:s"   => \$out_dir,
    "d|dir:s"      => \$dir,
    "f|function:s" => \$func_file,
    "n|name:s"     => \$output_id,
    "g|graphic:s"  => \$graphic_config_file,
    "help|h|?"     => sub { pod2usage( { -verbose => 2, -exitval => 1 } ) }
) || die("Cannot recognize usage. Please use --help to see options");

if ( !$file_in1 ) {

    pod2usage(2);
    exit();

}

#Initializes & runs the
main();

#This is the main sub reutine that performs all the task
sub main {

    if ($verb) {

        warn "Welcome to PanNav Creation Software\n\n";
        warn "The PanGenome File located in $file_in1 will be used..\n\n\n";

    }

    #Makes the output directories
    make_output_dirs();

    # This is the hash of the header of the tables
    #The default gene and region tabs
    $head{Region} = [ "ID", "Type", "# of Genes", "Associated Terms" ];
    $head{Gene} = [ "ID", "Name", "# of Genomes", "Region", "Associated Terms" ];

    #Setting the column widths of the file
    $col_width{Region} = [ .20, .20, .10, .50 ];
    $col_width{Gene} = [ .10, .20, .10, .20, .40 ];

    foreach my $a ( keys(%head) ) {

        for ( my $i = 0 ; $i < scalar( @{ $head{$a} } ) ; $i++ ) {

            $jsHead{$a}->{ $head{$a}->[$i] } = 1;

        }

    }

    #Checks to see if a tree is included, and if so loads it in
    #If the tree is loaded in, then tree_level is set to the defaut value (4): when it is not zero this means
    #that a tree is included
    if ($tree_file) {

        warn "Using $tree_file to incorporate tree in visualization\n\n";
        make_tree_json($tree_file);

    }

    #This is the default function if and only if there is a function file to will assign the rest of the
    $defaultFunct = "Hypothetical";
    if ( !$func_file ) { $defaultFunct = "NA"; }    #This is the case where no functions have been assigned. All will be NA
    $func_list{$defaultFunct} = 1;                  #This is the default assigment

    #file list;
    if ( -e $file_in1 ) {

        if ($verb) { warn "Reading in PanGenome information...\n"; }
        get_gene_info_new();

        if ($verb) { warn "Finished\n\n\n"; }

    } else {

        die("Cannot find $file_in1 for reading in files. Please double check...\n\n");

    }

    $rgb{$defaultFunct} = "dddddd";    #Default Function color is grey
    $rgb{"Background"}  = "e0e0e0";    #Background color is the
    $rgb{"fGR"}         = "606060";    #Default fGR color is red: should I change this?

    #Reading in functional information and annotation
    if ($func_file) {

        if ($verb) {

            warn "Reading in gene function information...\n";

        }
        read_function_file();
        #
        foreach my $funcStr (@funcInfo) {

            load_functions($funcStr);

        }
        if ($verb) {

            warn "Finished ...\n\n\n";

        }

    }

    #If the user wants to change some of the graphic size and color values
    if ($graphic_config_file) {

        read_config_file();

    }

    if ($verb) {

        warn "Starting the main image..\n";

    }

    #Write the javascript files used by the scripts: main (used by home page), fgi (fgi/core pages), and gene page
    fgi_javascript();
    make_main_javascript();
    write_geneFile_javascript();

    #making the main figure: the sub-script also uses this secondary and tertiary pages
    make_main_figure();
    if ($verb) {

        warn "Finished all\n";

    }

}

#All the files used in the image

sub make_main_figure() {

    #This is the main html file. Already the directory has been created
    open( SVG_MAIN, ">", $out_dir . "/$output_id/main.html" ) or die();

#%out is the main variable to store all the SVG information for the main pages images with the chromosome number
#as the hash index.
#Since SVGs doen't have z scores, we use staggered hash format for each chromosome with 100 being the highest, ie first drawn
#
#Number 100 is the maximum level allowed
#Number 100 is the back ground
#Number 2 is for the preview windows. It is also in a seperate hash key "preview"
#Number 1 is for the central ring + region
#Only levels 100 and 1 are used in the thumbnail

    #Other keys include:
    #key "load" shows the loading screen (currenly rotating exapanding and diminishing JCVI square logo)
    #key "script" loads the javascript
    #There are two scripts: first which loads the default variables and second which calls the script
    #Key "beg" is the opening for the HTML: this includes the javascript, legend, save
    #Key "end" is the closing for the HTML
    #Key "load" is the animation for the loading screen
    #Key "script" is the script
    #Key "table" is the table
    #Key "thumb" contains the thumbnail sketches

    my %out;

    #The ratio of the thumbnail to the main images
    my $assem_core_num = scalar( keys(%cores) );
    my %thumb_list;

    #Making the html images. Ordered similarly as %out -->
    my %html;

    $html{beg} = "<!DOCTYPE html><html lang=\"en\">";

    #getting the javascript functions
    $html{beg} .= "<head> <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"/>";

    #Initializeing the image by drawing
    %out = &svg_init_image( 2 * ( $border + $max_radius ), 2 * ( $border + $max_radius ), \%out );
    #
    $out{beg} = sprintf(
"<div id=\"mainDiv\" style=\"visibility: hidden; left:%fpx; top:0px; position:absolute;\" onmousemove=\"mouseMove()\">",
        1 * $gene_height )
      . $out{beg};

    #This is the loading screen described above. The blocks are the squares and twirl is the central screen
    $out{load} = make_loading_screen();

    #Setting the loading of the initial javascript variables
    my $background = $rgb{"Background"};    #setting the background color for the arc

    my $core_num = 1;                       #starting the number of the core genes

    #This is making the save disk icons for SVG and PNGs
    $out{beg} .= sprintf(
        "<path d=\"M%f %f L%f %f A%f %f 0 0 1 %f %f Z\" fill=\"lightgray\"/>\n",
        ( $border + $max_radius ),
        ( $border + $max_radius ),
        $border * 2,
        ( $border + $max_radius ),
        ( $max_radius - $border ),
        ( $max_radius - $border ),
        2 * $max_radius,
        $border + $max_radius
    );

 #Making the help button that takes you to the help page. Currently this is the manual PDF on the github site
 $out{beg} .= sprintf("<a xlink:title=\"Help\" target=\"_blank\" href=\"https://github.com/JCVenterInstitute/PanACEA/blob/master/PanACEA.manual.pdf\"><svg id = \"helpButton\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\"><path d=\"M29.898 26.5722l-4.3921 0c-0.0118,-0.635 -0.0177,-1.0172 -0.0177,-1.1583 0,-1.4229 0.2352,-2.5929 0.7056,-3.5102 0.4704,-0.9231 1.417,-1.952 2.8281,-3.1044 1.4111,-1.1465 2.2578,-1.8991 2.5282,-2.2578 0.4292,-0.5585 0.6409,-1.1818 0.6409,-1.8579 0,-0.9408 -0.3763,-1.7463 -1.1289,-2.4224 -0.7526,-0.6703 -1.7639,-1.0054 -3.0397,-1.0054 -1.2289,0 -2.2578,0.3527 -3.0868,1.0524 -0.8232,0.6997 -1.3935,1.7698 -1.7051,3.2044l-4.4391 -0.5527c0.1234,-2.0578 0.9995,-3.8041 2.6223,-5.2387 1.6286,-1.4346 3.757,-2.152 6.4029,-2.152 2.7752,0 4.9859,0.7291 6.6322,2.1814 1.6404,1.4522 2.4635,3.1397 2.4635,5.0741 0,1.0642 -0.3057,2.0755 -0.9054,3.028 -0.6056,0.9525 -1.8933,2.2519 -3.8688,3.8923 -1.0231,0.8525 -1.6581,1.5346 -1.905,2.052 -0.2469,0.5174 -0.3587,1.4405 -0.3351,2.7752zm-4.3921 6.5087l0 -4.8389 4.8389 0 0 4.8389 -4.8389 0z\" style=\"stroke: black; fill: yellow;\"></svg></a>",
	( $border + $max_radius - 80 ), 
	( $border + $max_radius - $border * 1.90 -5),
	40,
	40);
 #Making the Disk image by: (1) shrinking it with $trans; (2) adding the action; and (3) making a new variable containing the
 #This gives the location(translate) and appropriate size for the main page save SVG image
    my $trans = sprintf(
        "translate(%f,%f) scale(%f, %f)",
        ( $border + $max_radius - 20 ),
        ( $border + $max_radius - $border * 1.90 ),
        0.045, 0.045
    );
    my $action     = "onclick=\"saveSVG(\'svg\')\"";
    my $disk_image = $disk_image_orig;
    $disk_image =~ s/TRANS/$trans/g;
    $disk_image =~ s/ACTION/$action/g;
    $disk_image =~ s/TEXT/SVG/g;
    $disk_image =~ s/FS/250/g;

    #Adding disk image to the SVG
    $out{beg} .= $disk_image;

    #Same as above but for PNG
    $trans = sprintf(
        "translate(%f,%f) scale(%f, %f)",
        ( $border + $max_radius + 20 ),
        ( $border + $max_radius - $border * 1.90 ),
        0.045, 0.045
    );
    $action     = "onclick=\"saveSVG(\'png\')\"";
    $disk_image = $disk_image_orig;
    $disk_image =~ s/TRANS/$trans/g;
    $disk_image =~ s/ACTION/$action/g;
    $disk_image =~ s/TEXT/PNG/g;
    $disk_image =~ s/FS/250/g;

    $out{beg} .= $disk_image;

    #Legend Title
    $out{beg} .= "<g id=\"Legend1\" visiblity=\"visible\">\n";
    $out{beg} .= sprintf(
"<text id=\"LegendText\" x=\"%f\" y=\"%f\" fill=\"black\" font-size=\"40\" text-anchor=\"middle\" numChrom=\"$core_cnt\" onclick=\"makeThumbnails(evt)\">%s</text>\n",
        ( $border + $max_radius ),
        ( $border + $max_radius - $border * 1.95 ), "Legend"
    );

    #Legend: The different functions
    my @curr_y;    #Keeps track of current y position of the function in the table
    my @curr_x;    #Keeps track of current x position of the function in the table

    my $num_leg_cols = 3;    #Setting the number of columns in the legend

    #Making the starting x and y position for each column
    my $x1 = ( ( $max_radius - $border ) * sin( deg2rad(310) ) ) + ( $border + $max_radius );
    my $y1 = ( ( $max_radius - $border ) * cos( deg2rad(310) ) * -1 ) + ( $border + $max_radius );
    for ( my $i = 0 ; $i < $num_leg_cols ; $i++ ) {

        $curr_x[$i] = $x1 + ( $i * $border * 1.15 );
        $curr_y[$i] = $y1;

    }
    my $i = 0;
    if ($verb) {

        warn "Found ", scalar( keys(%func_list) ), " different functions\n\n";

    }

    #func_list: this cycles through by column, and once it
    foreach my $a ( keys(%func_list) ) {

        $x1 = $curr_x[$i];
        $y1 = $curr_y[$i];

        #in case of long names, this function (1) switchs words for equivilant symbols (ie and -> &) and making multi-line
        my $id_name = split_id( $a, $x1 + $gene_height * 0.6, 22 );

        #Getting the number of lines
        my $n = -1;
        while ( $id_name =~ /\<\/tspan\>/g ) { $n++; }

        #Adding the name, the color rectangle. All should be clickable to turn it on/off
        $out{beg} .= sprintf( "<g onclick=\"runType(\'GeneType\',evt,\'\')\" href=\"%s\" fillcol=\"%s\">", $a, "#" . $rgb{$a} );

#The rectangle around: thought that I might add a box if it's clicked. Currently turrned off (instead all others fade when one is clicked)
        $out{beg} .= sprintf(
"<rect href=\"%s\" id=\"%s\" class=\"LegRect\" x = \"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"none\" visiblility=\"hidden\"/>",
            $a, $a . "rect",
            $x1, $y1, $gene_height * 0.75 + $border * 1.15, $gene_height
        );

        #The colored rectange before
        $out{beg} .= sprintf(
            "<rect  href=\"%s\" id=\"%s\" class=\"LegCol\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" fillcol=\"%s\"/>\n",
            $a, $a . "col", $x1,
            $y1 + $gene_height * ($n) / 5,
            $gene_height * 0.5,
            $gene_height / 3.0,
            "#" . $rgb{$a}, "#" . $rgb{$a}
        );

        #The text printed out. Can be multiple lines thanks to col_names
        $out{beg} .= sprintf(
"<text href=\"%s\" id=\"%s\" class=\"LegType\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"auto\" font-size=\"10\" dominant-baseline=\"hanging\" text-anchor=\"start\" fill=\"black\" fillcol=\"%s\">%s</text>\n",
            $a,
            $a . "Text",
            $x1 + $gene_height * 0.6,
            $y1, $border * 1.15, "#" . $rgb{$a}, $id_name
        );
        $out{beg} .= sprintf("</g>");

        $y1 += $gene_height * ( ( $n + 1 ) / 2.5 + 0.15 );    #cycling through columns
        $curr_y[$i] = $y1;
        $i++;
        if ( $i == 3 ) { $i = 0; }

    }

    $out{beg} .= "</g>";

    ####################################################################################
    # Cycling through the pan-chromosome by:
    #		1. Looking at each core gene and:
    #			1.1 If it is a core gene,
    #				1.1.1 add information to cor_reg hash and goto 1
    #			1.2 If it it is a puedo-core gene (i.e. flexible region)
    #				1.2.1 If core genes in the cor_reg hash:
    #					1.2.1.1 Make core region preview pane and page from cor_reg
    #					1.2.1.2 Draw core region clickable segment on the chromosomes
    #					1.2.1.2 Reset cor_reg hash
    #				1.2.2 Find all genes in associated with fGRs and build fgi_reg hash
    #				1.2.3 Make fGR region preview pane and page from fgi_reg
    #		2. If core genes in cor_reg hash
    #			2.2.1 Make core region preview pane and page from cor_reg
    #			2.2.2 Draw core region clickable segment on the chromosomes
    ####################################################################################

    my %cor_reg;    #Hash with a list of core region information passed to the svg draw page
                    #Keys: 	st -> starting bp
                    #		end -> ending bp
                    #		coords -> array with gene

    my %reg2gene;   #Hash that maps the region to a gene

    my $core_type     = "";    # contains the annotations for each core gene: used for table
    my $core_oth      = "";    # contains any additional annotation for each core gene; ibid
    my $cur_gene_list = "";    #Contains a semi-colon sperated list of gene IDs for each region
    my $old_st        = 0;     #Used to keep track of the position of the core region start
    my $even          = 0;     # used to iterate the position of the core regions to shuffle the verticle postion of the fGRs

    for my $cur_chr ( keys(%cores) ) {

        $out{script_chr} .= sprintf( "\n\nvar %s = \'{", "fastaJSONPlasmid" . $cur_chr );

        if ($verb) {

            warn "Starting the image for Core $cur_chr..\n";

        }

        my $c3 = 0;            #Number of chromosomes in the thumbnails
        if ( $thumb_list{Assembly_core} ) { $c3 = scalar( @{ $thumb_list{Assembly_core} } ); }
        $thumb_list{Assembly_core}->[$c3] = $cur_chr;

        if ( $cores{$cur_chr}->{type} eq "cycle" ) {

            #Make a complete circle
            $cur_st = 0;

            #Setting the size of the chromosome. This is to make sure that the ticks/genes are correctly spaced
            $seq_len = $cores{$cur_chr}->{sz};

            #Drawing the ticks around the chromosome
            #Trying to get ticks that have have intervals of a factor of 10 and number between 10 to 100
            my $tick_space = 10**( int( log($seq_len) / log(10) ) - 1 );
            $out{$cur_chr} = &place_svg_ticks( $tick_space, $seq_len, $out{$cur_chr} );

            #Drawing the background of a full circle
            $out{$cur_chr} = &svg_draw_arc_full(
                100, &convert_layer(1),
                &convert_layer(1) + $gene_height,
                "#" . $background,
                0, 0, $out{$cur_chr}
            );

        }
        if ( $cores{$cur_chr}->{type} eq "chain" ) {

            #Change to incomplete circle- 300 degree in a horseshoe shape
            $cur_st = $cores{$cur_chr}->{sz} * 7.5 / 300;

            #This is to make sure that the ticks/genes are in the appropiate place
            #The functions are designed to work on a full
            $seq_len = $cores{$cur_chr}->{sz} * 315 / 300;

            #svg_draw_arc_seg-> change starting position on arc but not genome.
            $out{$cur_chr} = &svg_draw_arc_seg(
                100, 0, $cores{$cur_chr}->{sz},
                &convert_layer(1),
                &convert_layer(1) + $gene_height,
                "#" . $background,
                "", 0, "", "", $out{$cur_chr}
            );

            #Trying to get ticks that have have intervals of a factor of 10 and number between 10 to 100
            my $tick_space = 10**( int( log($seq_len) / log(10) ) - 1 );

            #A modified verion of the above
            $out{$cur_chr} = &place_svg_ticks_arc( $tick_space, $cores{$cur_chr}->{sz}, $out{$cur_chr} );

        }
        my $terms;
        my %term_cnt;
        my $out_cnt = 0;

        #Counting through the core genes: can also include psuedo-core genes, that is fGRs
        foreach my $core ( @{ $core_list{$cur_chr} } ) {

#Each corde gene includes the chromosome ($chr), the ID, the starting BP, the ending BP, the long gene name, the type (core/fgi), and the length in BP
            my ( $chr, $ID, $st, $end, $def, $type, $len ) = @$core;

            #Setting the start value for the next group
            $old_st = $end + 1;

            if ( $type eq "fGI" ) {

                #Getting the heights of the Regions on the pan-chromosome
                my $h1 = &convert_layer(1);
                my $h2 = $h1 + ( $gene_height * 0.75 );

                #This gets the fGR staggering variable
                if ($even) { $h1 += $gene_height * .25; $h2 += $gene_height * .25; $even = 0; }
                else       { $even = 1; }

                #Checking to see if there is a core-region (ie area with core-genes)
                if (%cor_reg) {

                    #The location of the Core Region Detail Page
                    my $new = "/core/Core_Region.$core_cnt.html";
                    my $c   = 0;                                    #The row of the Region in the Region Table Array
                    if ( $jso{Region} ) { $c = scalar( @{ $jso{Region} } ); }

                    #The values to be added to the Region Table
                    my @tmp = (
                        "Core_Region.$core_cnt", "Core Region", scalar( @{ $cor_reg{cords} } ),
                        $core_type, "CORE" . $core_num
                    );
                    for ( my $i = 0 ; $i < scalar( @{ $head{Region} } ) ; $i++ ) {

                        $jso{Region}->[$c]->{ $head{Region}->[$i] } = $tmp[$i];

                    }
                    $reg2gene{ "CORE" . $core_num } = $cur_gene_list;

                    #Adding information to json array for the Region to be used in the javascript but not shown in the tba
                    $jso{Region}->[$c]->{type_ref}   = $core_type;
                    $jso{Region}->[$c]->{oth_ref}    = $core_oth;
                    $jso{Region}->[$c]->{core_clust} = $cur_chr;
                    $jso{Region}->[$c]->{href}       = "CORE" . $core_num;

                    #Counting the number of regions as so to name them
                    $core_cnt++;

                    #Drawing the core region preview pane
                    #$tmp is a string holding the SVG and then the HTML values for the pane and then the page
                    my $tmp;
                    $tmp = &make_core_region_svg( $new, "core" . $core_cnt, \%cor_reg, \%ont, $tmp );

                    #Adding the preview pane to the SVG $out variable in the preview key at level 2
                    $out{$chr}->{preview}->{2} .= $tmp;

                    #Making the core region detail page
                    $tmp = &make_core_region_page( "core" . $core_cnt, \%cor_reg, $tmp );
                    $tmp =~ s/hidden/visible/g;    #Probably can change this
                    open( TMP_OUT, ">", $out_dir . "/$output_id/" . $new );
                    print TMP_OUT $tmp;
                    close(TMP_OUT);

                    #Drawing the clickable area on the pan-chromosome
                    $out{$chr} = &svg_draw_arc_seg(
                        1, $cr_start, $st - 1, &convert_layer(1), &convert_layer(1) + $gene_height, "#" . $background, $new,
                        sprintf( "CORE REGION $cr_start to %d\n%d Core Genes", ( $st - 1 ), $cr_cnt ), "CORE" . $core_num,
                        "Region", $out{$chr}
                    );

                    #Counting the number of core regions
                    $core_num++;

                }

                #Resetting the core_region values
                $cur_gene_list = "";
                $core_type     = "";
                $core_oth      = "";

                #String with the file id of the
                my $new = $out_dir . "/$output_id/FGI.$core_cnt.html";

                #The FGI_region hash used to make the preview pane and the detail page
                #Keys: 	st -> starting bp
                #		end -> ending bp
                #		number -> count of the different gene presence and absences paths observed in the fGR
                #		order -> array with the different gene paths sorted by the number of genomes with that order
                #		cnt -> hash with key of fgi ID with number of genomes containing the path
                # 		The ID for most of these is the path ID
                #		st_id -> hash with key of fgi ID with the bounding core gene ID
                #		st_sz -> hash with key of fgi ID with the bounding core gene size
                #		st_dir -> hash with key of fgi ID with the bounding core gene direction
                #		st_col -> hash with key of fgi ID with the bounding core gene color
                #		end_id -> hash with key of fgi ID with the bounding core gene ID
                #		end_sz -> hash with key of fgi ID with the bounding core gene size
                #		end_dir -> hash with key of fgi ID with the bounding core gene direction
                #		end_col -> hash with key of fgi ID with the bounding core gene color
                #		gen -> string with semi-colon separated list of genomes with
                #		coords -> array with gene
                my %fgi_reg;

                my %fgi_ids;    #Keeps track of the gene IDs in each gene path using a simple 1 if its in the path
                $fgi_reg{st} = min( $st, $end );    #Ignoring the direction
                $fgi_reg{end} = max( $st, $end );
                $fgi_reg{number} = scalar( keys( %{ $list{$ID} } ) );    #Getting the counts of the genomes

             #@sort is the based on the number of genes in the path this is subsequently used to sort the genes on the x-axis
             #@sort2 is based on the number of genomes containing this path. This is the order used on the pane/page y-axis
                my @sort =
                  sort { scalar( @{ $list{$ID}->{$b}->{arr} } ) <=> scalar( @{ $list{$ID}->{$a}->{arr} } ) }
                  keys( %{ $list{$ID} } );
                my @sort2 = sort { $list{$ID}->{$b}->{cnt} <=> $list{$ID}->{$a}->{cnt} } keys( %{ $list{$ID} } );

                $fgi_reg{order} = \@sort2;

                #Going through each of the paths
                foreach my $a (@sort2) {

                    $fgi_reg{cnt}->{$a} = ( $list{$ID}->{$a}->{cnt} + 0 );

                    #Adding information about the starting bounding core gene if it exists for this gene path
                    #See list above for the variables being collected
                    if ( $list{$ID}->{$a}->{st} ) {

                        $fgi_reg{st_id}->{$a}  = $list{$ID}->{$a}->{st};
                        $fgi_reg{st_sz}->{$a}  = abs( $core_size{ $list{$ID}->{$a}->{st} } );
                        $fgi_reg{st_dir}->{$a} = 1;

                        if ( $core_size{ $list{$ID}->{$a}->{st} } < 0 ) {

                            $fgi_reg{st_dir}->{$a} = -1;

                        }
                        if ( $gene_func{ $list{$ID}->{$a}->{st} } ) {

                            $fgi_reg{st_col}->{$a} = $rgb{ $gene_func{ $list{$ID}->{$a}->{st} } };

                        } else {

                            #Setting empty/blank color values to the default function
                            $fgi_reg{st_col}->{$a} = $rgb{$defaultFunct};

                        }

                    } else {

                        #The blank start core gene/break values
                        $fgi_reg{st_id}->{$a}  = 0;
                        $fgi_reg{st_sz}->{$a}  = 0;
                        $fgi_reg{st_dir}->{$a} = 0;
                        $fgi_reg{st_col}->{$a} = "ffffff";

                    }
                    if ( $list{$ID}->{$a}->{end} ) {

                        #Adding information about the starting bounding core gene if it exists for this gene path
                        #See list above for the variables being collected

                        $fgi_reg{end_id}->{$a}  = $list{$ID}->{$a}->{end};
                        $fgi_reg{end_sz}->{$a}  = abs( $core_size{ $list{$ID}->{$a}->{end} } );
                        $fgi_reg{end_dir}->{$a} = 1;
                        if ( $gene_func{ $list{$ID}->{$a}->{end} } ) {

                            $fgi_reg{end_col}->{$a} = $rgb{ $gene_func{ $list{$ID}->{$a}->{end} } };

                        } else {

                            $fgi_reg{end_col}->{$a} = $rgb{$defaultFunct};

                        }

                    } else {

                        #No bounding end core gene
                        $fgi_reg{end_id}->{$a}  = 0;
                        $fgi_reg{end_sz}->{$a}  = 0;
                        $fgi_reg{end_dir}->{$a} = 0;
                        $fgi_reg{end_col}->{$a} = "ffffff";

                    }
                    $fgi_reg{end_type}->{$a} = $gene_func{ $list{$ID}->{$a}->{end} };
                    $fgi_reg{st_type}->{$a}  = $gene_func{ $list{$ID}->{$a}->{st} };
                    $fgi_reg{gen}->{$a}      = $list{$ID}->{$a}->{gen};

                }

                my %mat;

                #mat is hash to store the array with the genes so it can be ordered and displayed
                #keys are arr -> array with the genes in order; loc -> hash where the gene id points to its order

                #For the first gene in the fgi

                for ( my $j = 0 ; $j < scalar( @{ $list{$ID}->{ $sort[0] }->{arr} } ) ; $j++ ) {

                    $mat{arr}->[$j] = $list{$ID}->{ $sort[0] }->{arr}->[$j];
                    my $a = $list{$ID}->{ $sort[0] }->{arr}->[$j];
                    if ($a) {

                        if ( !$clusters{$a} ) {

                            die("Cannot find a cluster for $a. Ending PanACEA...\n\n");

                        }

#Adding the region information to the table. This could allow for clicking on a gene to highlight it's location in the region
                        if (! $jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{$head{Gene}->[3]})
						{
							$jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{ $head{Gene}->[3] } = $ID;
                        }
						else
						{
							$jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{ $head{Gene}->[3] } .= ",". $ID;
						}
						$jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{core_clust}         = $cur_chr;
                        $jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{href}               = "CL_" . $list{$ID}->{ $sort[0] }->{arr}->[$j];
						$jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{detail} =
                          "CL_" . $list{$ID}->{ $sort[0] }->{arr}->[$j];

                        $cur_gene_list .= "CL_" . $list{$ID}->{ $sort[0] }->{arr}->[$j] . ";";
                        $mat{loc}->{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } = $j + 1;

                        #Adding values to the functional annotation of the Region, both the type ref and the other ref
                        while (
                            $jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{type_ref} =~ /([^\;]+)/g ) {

                            my $in1 = $1;
                            if ( $core_type !~ /$in1/ ) { $core_type .= $in1 . ";"; }

                        }
                        while (
                            $jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{oth_ref} =~ /([^\;]+)/g ) {

                            my $in1 = $1 . ";";
                            if ( $ont{$in1} ) {

                                $term_cnt{$in1}++;
                                $terms->{$in1} = $ont{$in1};

                            }
                            if ( !$term2type{$in1} ) {

                                my $type_ref =
                                  $jso{Gene}->[ $gene_num{ $list{$ID}->{ $sort[0] }->{arr}->[$j] } ]->{type_ref};
                                $term2type{$in1} = $type_ref;
                                if ( !$type2term{$type_ref} ) {

                                    $type2term{$type_ref}->[0] = $in1;

                                } else {

                                    push @{ $type2term{$type_ref} }, $in1;

                                }

                            }
                            if ( $core_oth !~ /$in1/ ) { $core_oth .= $in1; }

                        }

                        #Flags the gene as being in the fgi
                        $fgi_ids{ $list{$ID}->{ $sort[0] }->{arr}->[$j] }->{ $sort[0] } = 1;

                    }

                }

                #For all subsequent genes in the fGI
                for ( my $i = 1 ; $i < scalar(@sort) ; $i++ ) {

                    #Prev is the order value, which is used to insert a previously unseen gene into the order
                    my $prev = 0;

                    #For each of the paths with the Gene
                    for ( my $j = 0 ; $j < scalar( @{ $list{$ID}->{ $sort[$i] }->{arr} } ) ; $j++ ) {

                        my $a = $list{$ID}->{ $sort[$i] }->{arr}->[$j];
                        if ($a) {

                            $fgi_ids{$a}->{ $sort[$i] } = 1;

                            #If there is an order position for the gene, use it
                            if ( $mat{loc}->{$a} ) {

                                $prev = $mat{loc}->{$a};

                            } else {

                                #Otherwise, cycle through all the previously seen genes to get a count
                                for ( my $k = scalar( @{ $mat{arr} } ) ; $k > ($prev) ; $k-- ) {

                                    $mat{arr}->[$k] = $mat{arr}->[ $k - 1 ];
                                    $mat{loc}->{ $mat{arr}->[$k] }++;

                                }
                                $mat{arr}->[$prev] = $a;
                                $mat{loc}->{$a} = $prev + 1;
                                $prev++;

                            }

                        }

                    }

                }

                #If the matrix array exists
				my $fgr_gene_cnt = 0;
                
                if ( $mat{arr} ) {

					$fgr_gene_cnt = scalar( @{ $mat{arr} } );               
				   #Go through the array of genes in the fgi
                    for ( my $i = 0 ; $i < scalar( @{ $mat{arr} } ) ; $i++ ) {

                        my $a = $mat{arr}->[$i];
                        if ($a) {

                            $fgi_grps{ $fgi_member{$a} }->{hit} = $ID;
                            if ( !$clusters{$a} ) { die("$ID $i $a"); }

                            #Put the gene information into the coordinate array
                            push(
                                @{ $fgi_reg{cords} },
                                [
                                    $i,                                                $clusters{$a}->{num_of_members},
                                    $clusters{$a}->{protein_name},                     $fgi_member{$a},
                                    $fgi_ids{$a},                                      "CL_$a",
                                    $rgb{ $jso{Gene}->[ $gene_num{$a} ]->{type_ref} }, $clusters{$a}->{start},
                                    $clusters{$a}->{end},                              $jso{Gene}->[ $gene_num{$a} ]
                                ]
                            );

                        }

                    }

                    my $c = 0;

                    #If there are any genes in the fGR...
                    if ( scalar( @{ $fgi_reg{cords} } ) > 0 ) {

                        #if there are- make the pane and the page
                        my $fig;

                        #drawing the fgi page and writing it to a html file
                        $fig = &make_fgi_page_svg( "fgi" . $ID, \%fgi_reg, \%ont, $fig );
                        $new = "fgi/$ID.html";

                        open( OUT, ">", $out_dir . "/$output_id/" . $new );
                        print OUT "$fig\n";
                        close(OUT);

                        #drawing the preview pane
                        $out{$chr}->{preview}->{2} .= &make_fgi_region_svg( $new, \%fgi_reg, \%ont, $out{$chr}->{2} );

                    }

                    #Counting the number of regions to keep in the table
                    if ( $jso{Region} ) { $c = scalar( @{ $jso{Region} } ); }

                    #Putting the region information into a table row for the table
                    my @tmp = ( "FGR.$core_cnt", "fGR", scalar( @{ $fgi_reg{cords} } ), $core_type, $ID );
                    for ( my $i = 0 ; $i < scalar( @{ $head{Region} } ) ; $i++ ) {

                        $jso{Region}->[$c]->{ $head{Region}->[$i] } = $tmp[$i];

                    }

                    #More table information for put in the table jso
                    $jso{Region}->[$c]->{href}       = $ID;
                    $jso{Region}->[$c]->{core_clust} = $cur_chr;
                    $jso{Region}->[$c]->{oth_ref}    = $core_oth;

                    $jso{Region}->[$c]->{type_ref} = $core_type;
                    $reg2gene{$ID} = $cur_gene_list;

                    #Drawing the fGI on the pan-chromosome map
                    $core_type = "";
                    $core_oth  = "";
                    $out{$chr} =
                      &svg_draw_arc_seg( 1, $st, $end, $h1, $h2, "#" . $rgb{fGR}, $new, sprintf( "FGR REGION %d\n%d Total Genes", ( $core_cnt ), $fgr_gene_cnt ), $ID, "Region", $out{$chr} );

                }

                #incrementing the core region count and resetting
                $cr_start = $end + 1;
                $cr_cnt   = 0;
                %cor_reg  = ();

            } elsif ( $type eq "CL" ) {

                #for the Core genes, just add the gene information to cor_reg

                #add to the number of genes in the core region
                $cr_cnt++;

                #Getting the starting Base pair
                if ( !%cor_reg ) {

                    $cor_reg{st} = &min( $st, $end );

                }

                #Get the number from the Core Gene ID
                $ID =~ /_(\d+)/;
                my $num1 = $1;

                #Adding id to the gene list for the core region
                $cur_gene_list .= "$ID;";

                #Adding information the core gene infromation to the table JSON
                $jso{Gene}->[ $gene_num{$num1} ]->{href}               = $ID;                  #"CORE" . $core_num;
                $jso{Gene}->[ $gene_num{$num1} ]->{detail}             = $ID;
                $jso{Gene}->[ $gene_num{$num1} ]->{core_clust}         = $cur_chr;
                $jso{Gene}->[ $gene_num{$num1} ]->{ $head{Gene}->[3] } = "CORE" . $core_num;

                #adding all the TYPE ids of the gene to the core type
                while ( $jso{Gene}->[ $gene_num{$num1} ]->{type_ref} =~ /([^\;]+)/g ) {

                    my $in1 = $1;
                    if ( $core_type !~ /$in1/ ) { $core_type .= "$in1;"; }

                }

                #adding all the reference (such as GO and ARO Terms) ids of the gene to the core reference type
                while ( $jso{Gene}->[ $gene_num{$num1} ]->{oth_ref} =~ /([^\;]+)/g ) {

                    my $in1 = $1 . ";";
                    if ( $core_oth !~ /$in1/ ) { $core_oth .= "$in1"; }

                }

                push @$core, "#" . $rgb{ $jso{Gene}->[ $gene_num{$num1} ]->{type_ref} };
                push @$core, $jso{Gene}->[ $gene_num{$num1} ];

                push( @{ $cor_reg{cords} }, $core );

                $cor_reg{end} = max( $st, $end );
                my $same = $ID . "detail";

                #Also draw the appropriate line on the inner two rings

                my $genomeList    = "";    #string list of the genomes which contains a core gene
                my $genomeLineLen = 0;     #x-location of the genome list. Used to add breaks to the genome list
                foreach my $genomeID ( keys( %{ $seqs{byClust}->{$num1} } ) ) {

                    $genomeLineLen += length($genomeID);
                    if ( $genomeLineLen > 80 ) {

                        $genomeList .= "<br>";
                        $genomeLineLen = length($genomeID);

                    }                      #Adding the break if needed in the line is > 80 characters
                    $genomeList .= "$genomeID,";

                }

                $out{script_chr} .= sprintf(
"\"%s\":{\"Cluster\":\"%s\", \"Name\":\"%s\",\"# of Genomes\":\"%s\",\"Functional IDs\":\"%s\",\"Associated Terms\":\"%s\",\"Sequence\":\"%s\",\"Maximum Length\":\"%f\", \"Minimum Length\":\"%f\", \"Mean\":\"%f\", \"Standard Deviation\":\"%f\", \"Genomes\":\"%s\"},",
                    $ID,                                          $ID,
                    $clusters{$num1}->{protein_name},             $clusters{$num1}->{num_of_members},
                    $jso{Gene}->[ $gene_num{$num1} ]->{type_ref}, $jso{Gene}->[ $gene_num{$num1} ]->{oth_ref},
                    $seqs{cent}->{$num1}->{seq},                  $geneLenInfo{$num1}->{minLen},
                    $geneLenInfo{$num1}->{maxLen},                $geneLenInfo{$num1}->{mean},
                    $geneLenInfo{$num1}->{sd},                    $genomeList
                );

                my @term_list = split ";", $jso{Gene}->[ $gene_num{$num1} ]->{oth_ref};
                foreach my $a (@term_list) {

                    if ( $ont{$a} ) { $terms->{$a} = $ont{$a}; $term_cnt{$a}++; }

                }
                $out_cnt++;

                #Drawing the core gene line on the interior of the pan-chromosome circle
                if ( $st < $end ) {

                    #This is 5'
                    my $h1 = &convert_layer(2);
                    my $h2 = $h1 + ($gene_height);
                    $out{$chr} =
                      &svg_draw_arc_seg( 3, $st, $end, $h1, $h2, "#" . $rgb{ $jso{Gene}->[ $gene_num{$1} ]->{type_ref} },
                        $same, $ID, $cur_chr, "Gene", $out{$chr} );

                } else {

                    my $tmp = $st;
                    $st  = $end;
                    $end = $tmp;
                    my $h1 = &convert_layer(3);
                    my $h2 = $h1 + ($gene_height);
                    $out{$chr} =
                      &svg_draw_arc_seg( 3, $st, $end, $h1, $h2, "#" . $rgb{ $jso{Gene}->[ $gene_num{$1} ]->{type_ref} },
                        $same, $ID, $cur_chr, "Gene", $out{$chr} );

                }

            }

        }

        #If a core region is still not yet drawn at the end of the chromosome...
        if (%cor_reg) {

            #Draw it! This is the same as drawing it at beginning of an fGI
            my $chr = $cur_chr;
            my $new = "/core/Core_Region.$core_cnt.SVG";
            my $c   = 0;
            if ( $jso{Region} ) { $c = scalar( @{ $jso{Region} } ); }
            my @tmp =
              ( "Core_Region.$core_cnt", "Core Region", scalar( @{ $cor_reg{cords} } ), $core_type, "CORE" . $core_num );
            for ( my $i = 0 ; $i < scalar( @{ $head{Region} } ) ; $i++ ) {

                $jso{Region}->[$c]->{ $head{Region}->[$i] } = $tmp[$i];

            }

            $jso{Region}->[$c]->{type_ref}   = $core_type;
            $jso{Region}->[$c]->{oth_ref}    = $core_oth;
            $jso{Region}->[$c]->{core_clust} = $cur_chr;
            $jso{Region}->[$c]->{href}       = "CORE" . $core_num;
            $core_cnt++;

            my $tmp;
            $tmp = &make_core_region_svg( $new, 'core'.$core_cnt, \%cor_reg, \%ont, $tmp );
            $out{$chr}->{preview}->{2} .= $tmp;
            $tmp = &make_core_region_page( $new, \%cor_reg, $tmp );
            $tmp =~ s/hidden/visible/g;
            $reg2gene{ "CORE" . $core_num } = $cur_gene_list;

            open( TMP_OUT, ">", "$out_dir/$output_id/$new" );
            print TMP_OUT $tmp;
            close(TMP_OUT);
            $out{$chr} = &svg_draw_arc_seg(
                1, $cr_start, $old_st - 1, &convert_layer(1), &convert_layer(1) + $gene_height, "#" . $background, $new,
                sprintf( "CORE REGION $cr_start to %d\n%d Core Genes", ( $old_st - 1 ), $cr_cnt ), "CORE" . $core_num,
                "Region", $out{$chr}
            );
            $core_num++;
            $core_type     = "";
            $core_oth      = "";
            $cur_gene_list = "";

            #resetting the core
            %cor_reg = ();

        }
        chop( $out{script_chr} );
        $out{script_chr} .= "}\';\n";
        $out{script_chr} .= sprintf( "\n\nvar %s = \'{", "termJSONPlasmid" . $cur_chr );
        foreach my $a ( keys(%$terms) ) { $out{script_chr} .= sprintf( "\"%s\":\"%s\",", $a, $terms->{$a} ); }
        if ( $out_cnt > 0 ) {

            chop( $out{script_chr} );

        }
        $out{script_chr} .= "}\';\n";
        $out{script_chr} .= sprintf( "\n\nvar %s = \'{", "termCountJSONPlasmid" . $cur_chr );
        foreach my $a ( keys(%term_cnt) ) { $out{script_chr} .= sprintf( "\"%s\":\"%s\",", $a, $term_cnt{$a} ); }
        if ( $out_cnt > 0 ) {

            chop( $out{script_chr} );

        }

        $out{script_chr} .= "}\';\n";

    }

    #Creating a javacscript in the HTML file in the out hash script key
    $out{script} .= "<script type=\"text/javascript\" >
				//<![CDATA[";

    #Going through the fgis
    for my $cur_chr ( keys(%fgi_grps) ) {

        #If the fgi hasn't already been included
        if ( !$fgi_grps{$cur_chr}->{hit} ) {

            #If it is a circular fgi (plasmid?) write it to main svg
            if ( $fgi_grps{$cur_chr}->{type} eq "cycle" ) {

                my %cor_reg;
                my $c3 = 0;

                #Getting the number of Cycle fGRs to the thumbnail menu list
                if ( $thumb_list{"Cycle fGRs"} ) {

                    $c3 = scalar( @{ $thumb_list{"Cycle fGRs"} } );

                }
                $thumb_list{"Cycle fGRs"}->[$c3] = $cur_chr;

                #Adding all the fGR gene info region to the core list
                foreach my $core ( @{ $core_list{$cur_chr} } ) {

                    my ( $chr, $ID, $st, $end, $def, $type, $len ) = @$core;
                    if ( !%cor_reg ) {

                        $cor_reg{st} = min( $st, $end );

                    }
                    $cor_reg{end} = max( $st, $end );

                }
                my $st_dif = $cor_reg{st};
                $cor_reg{end} -= $cor_reg{st};
                $cor_reg{st} = 0;
                my $size = $cor_reg{end};    #size of the fGR
                $seq_len = $size;
                $cores{$cur_chr}->{sz} = $size;

                #drawing the circular region by drawing the arcs and the ticks
                $out{$cur_chr}->{beg} = sprintf(
                    "<svg height=\"%f\" width=\"%f\" id=\"Plasmid_%f\">\n",
                    2 * ( $border + $max_radius ),
                    2 * ( $border + $max_radius ), $cur_chr
                );
                $out{$cur_chr} = &place_svg_ticks( 10**( int( log($seq_len) / log(10) ) - 1 ), $seq_len, $out{$cur_chr} );
                $out{$cur_chr} = &svg_draw_arc_full(
                    100, &convert_layer(1),
                    &convert_layer(1) + $gene_height,
                    "#" . $background,
                    0, 0, $out{$cur_chr}
                );

                #Add the gene information for each fGR region via javascript
                $out{script} .= sprintf( "\n\nvar %s = \'{", "fastaJSONPlasmid" . $cur_chr );

                my $out_cnt;
                my $terms;
                my %term_cnt;
                foreach my $core ( @{ $core_list{$cur_chr} } ) {

                    my ( $chr, $ID, $st, $end, $def, $type, $len ) = @$core;

                    $ID =~ /_(\d+)/;
                    my $num1 = $1;
                    $cur_gene_list .= "$ID;";

                    #Adding information to the gene table JSON variable
                    $jso{Gene}->[ $gene_num{$num1} ]->{href}               = $ID;
                    $jso{Gene}->[ $gene_num{$num1} ]->{detail}             = $ID;
                    $jso{Gene}->[ $gene_num{$num1} ]->{core_clust}         = $cur_chr;
                    $jso{Gene}->[ $gene_num{$num1} ]->{ $head{Gene}->[3] } = "CORE" . $core_num;

                    #adding the gene type to the table
                    while ( $jso{Gene}->[ $gene_num{$num1} ]->{type_ref} =~ /([^\;]+)/g ) {

                        my $in1 = $1;
                        if ( $core_type !~ /$in1/ ) { $core_type .= "$in1;"; }

                    }

                    #adding the gene functional annotation to the table
                    while ( $jso{Gene}->[ $gene_num{$num1} ]->{oth_ref} =~ /([^\;]+)/g ) {

                        my $in1 = $1 . ";";
                        if ( $core_oth !~ /$in1/ ) { $core_oth .= "$in1"; }

                    }

                    push @$core, "#" . $rgb{ $jso{Gene}->[ $gene_num{$num1} ]->{type_ref} };
                    push( @{ $cor_reg{cords} }, $core );

                    my $h1 = &convert_layer(1);
                    my $h2 = $h1 + ($gene_height);

                    #Drawing the "plasmid" SVG:
                    $out{$cur_chr} = &svg_draw_plasmid_seg(
                        1,
                        $st - $st_dif,
                        $end - $st_dif,
                        $h1, $h2,
                        "#" . $rgb{ $jso{Gene}->[ $gene_num{$1} ]->{type_ref} },
                        "Plasmid" . $cur_chr,
                        $num1, $ID, "Gene", $out{$cur_chr}
                    );

                    #Adding gene information to the gene page information

                    my $genomeList    = "";    #string list of the genomes which contains a core gene
                    my $genomeLineLen = 0;     #x-location of the genome list. Used to add breaks to the genome list
                    foreach my $genomeID ( keys( %{ $seqs{byClust}->{$num1} } ) ) {

                        $genomeLineLen += length($genomeID);
                        if ( $genomeLineLen > 80 ) {

                            $genomeList .= "<br>";
                            $genomeLineLen = length($genomeID);

                        }                      #Adding the break if needed in the line is > 80 characters
                        $genomeList .= "$genomeID,";

                    }
                    $out{script} .= sprintf(
"\"%s\":{\"Cluster\":\"%s\", \"Name\":\"%s\",\"# of Genomes\":\"%s\",\"Functional IDs\":\"%s\",\"Associated Terms\":\"%s\",\"Sequence\":\"%s\"},\"Maximum Length\":\"%f\", \"Minimum Length\":\"%f\", \"Mean\":\"%f\", \"Standard Deviation\":\"%f\", \"Genomes\":\"%s\"",
                        $ID,                                          $ID,
                        $clusters{$num1}->{protein_name},             $clusters{$num1}->{num_of_members},
                        $jso{Gene}->[ $gene_num{$num1} ]->{type_ref}, $jso{Gene}->[ $gene_num{$num1} ]->{oth_ref},
                        $seqs{cent}->{$num1}->{seq},                  $geneLenInfo{$num1}->{minLen},
                        $geneLenInfo{$num1}->{maxLen},                $geneLenInfo{$num1}->{mean},
                        $geneLenInfo{$num1}->{sd},                    $genomeList
                    );

                    my @term_list = split ";", $jso{Gene}->[ $gene_num{$num1} ]->{oth_ref};
                    foreach my $a (@term_list) {

                        if ( $ont{$a} ) { $terms->{$a} = $ont{$a}; $term_cnt{$a}++; }

                    }
                    $out_cnt++;

                }

                $out{$cur_chr}->{end} .= "</svg>\n";
                chop( $out{script} );
                $out{script} .= "}\';\n";
                $out{script} .= sprintf( "\n\nvar %s = \'{", "termJSONPlasmid" . $cur_chr );
                foreach my $a ( keys(%$terms) ) { $out{script} .= sprintf( "\"%s\":\"%s\",", $a, $terms->{$a} ); }
                if ( $out_cnt > 0 ) {

                    chop( $out{script} );

                }
                $out{script} .= "}\';\n";
                $out{script} .= sprintf( "\n\nvar %s = \'{", "termCountJSONPlasmid" . $cur_chr );
                foreach my $a ( keys(%term_cnt) ) { $out{script} .= sprintf( "\"%s\":\"%s\",", $a, $term_cnt{$a} ); }
                if ( $out_cnt > 0 ) {

                    chop( $out{script} );

                }

                $out{script} .= "}\';\n";

            }

        }

    }
    $out{script} .= $out{script_chr};

    #After drawing the chromosomes, add the table information to the html as JSON strings able to be formatted

    #Drawing the table SVG
    #This is the blank bottom half of the circle that will be covered by the table
    $out{table} .= sprintf(
        "<path d=\"M%f %f L%f %f A%f %f 0 0 1 %f %f Z\" fill=\"#ffe4e1\"/>\n",
        ( $border + $max_radius ),
        ( $border + $max_radius ),
        2 * $max_radius,
        $border + $max_radius,
        ( $max_radius - $border ),
        ( $max_radius - $border ),
        $border * 2,
        ( $border + $max_radius )
    );

    my @table_types = keys(%jso);    #("Region", "Gene", "GO Term", "ARO Term");

    my $step = 110 / scalar(@table_types);
    my $st   = (305) / 180 * $PI;
    my $mid1 = ( 305 + 20 ) / 180 * $PI;
    my $mid2 = ( 305 + 90 ) / 180 * $PI;

    my $en = ( 305 + 110 ) / 180 * $PI;

    #light gray bottom part of the tabl
    $out{table} .= sprintf(
        "<path d=\"M%f %f A%f %f 0 0 0 %f %f L%f %f A%f %f 0 0 0 %f %f Z\" stroke=\"black\" fill=\"#bbffbb\" id=\"%s\" />\n",
        ( $border + $max_radius ) + sin($st) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($st) * ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($mid1) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid1) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($mid2) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid2) * ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($en) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($en) * ( $max_radius - $border * 1.5 ),
        "ButtonIDTableOn"
    );
    $out{table} .= sprintf(
"<text x=\"%f\" y=\"%f\" font-size=\"30\" text-anchor=\"middle\" alignment-baseline=\"middle\" id=\"tableButtonText\">Show Table</text></path>",
        $border + $max_radius,
        $border + $max_radius + ( $max_radius - $border * 2 )
    );
    $out{table} .= sprintf(
"<path d=\"M%f %f A%f %f 0 0 0 %f %f L%f %f A%f %f 0 0 0 %f %fZ\" stroke=\"black\" fill-opacity=\"0\" onclick=\"runType(\'TableStatus\',evt, \'\')\")/>",
        ( $border + $max_radius ) + sin($st) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($st) * ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($mid1) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid1) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($mid2) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid2) * ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($en) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($en) * ( $max_radius - $border * 1.5 )
    );

    #starting the table buttons
    $out{table} .= "<g id=\"tableButtons\" visibility=\"hidden\">\n";
    $out{table} .= sprintf(
        "<path d=\"M%f %f A%f %f 0 0 0 %f %f  Z\" stroke=\"#eeeeff\" fill=\"#eeeeff\" id=\"%s\" />\n",
        ( $border + $max_radius ) + sin($mid1) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid1) * ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($mid2) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid2) * ( $max_radius - $border * 1.5 ),
        ,
        "ButtonIDTableSave"
    );

    $trans = sprintf(
        "translate(%f,%f) scale(%f, %f)",
        $border + $max_radius,
        $border + $max_radius + ( $max_radius - $border * 1.75 ),
        0.045, 0.045
    );
    $action = "onclick=\"saveTableTxt()\"";

    #Adding disk image to the table
    $disk_image = $disk_image_orig;
    $disk_image =~ s/TRANS/$trans/g;
    $disk_image =~ s/ACTION/$action/g;
    $disk_image =~ s/TEXT/TXT/g;
    $disk_image =~ s/FS/250/g;
	my $saveTableDiskImage = $disk_image;
    $out{table} .= $disk_image;

#$out->{table} .= sprintf("<text x=\"%f\" y=\"%f\" font-size=\"20\" text-anchor=\"middle\" onclick=\"saveTableTxt()\" alignment-baseline=\"middle\" id=\"tableButtonText\">Save Table as TSV</text></path>", ;
    $out{table} .= sprintf(
        "<path d=\"M%f %f A%f %f 0 0 0 %f %f  Z\" stroke=\"#eeeeff\" fill=\"none\" id=\"%s\" onclick=\"saveTableTxt()\"/>\n",
        ( $border + $max_radius ) + sin($mid1) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid1) * ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + sin($mid2) * ( $max_radius - $border * 1.5 ),
        ( $border + $max_radius ) + cos($mid2) * ( $max_radius - $border * 1.5 ),
        ,
        "ButtonIDTableSave"
    );

    #Drawing all the type buttons on the table
    for ( my $i = 0 ; $i < scalar(@table_types) ; $i++ ) {

        $st = ( 305 + $step * $i ) / 180 * $PI;
        $en = ( 305 + $step * ( $i + 1 ) ) / 180 * $PI;
        my $id_in = $table_types[$i];
        $id_in =~ s/\A(\s)//g;    #warn "$id_in\n";

        #Add the svg portions incly
        $out{table} .= sprintf(
            "<defs><path id=\"%s\" d=\"M %f %f A%f,%f 0 0 0 %f,%f\"/></defs>\n",
            "textPath-" . $id_in,
            ( $border + $max_radius ) + sin($st) * ( $max_radius - $border * 1.25 ),
            ( $border + $max_radius ) + cos($st) * ( $max_radius - $border * 1.25 ),
            ( $max_radius - $border * 1.25 ),
            ( $max_radius - $border * 1.25 ),
            ( $border + $max_radius ) + sin($en) * ( $max_radius - $border * 1.25 ),
            ( $border + $max_radius ) + cos($en) * ( $max_radius - $border * 1.25 )
        );
        $out{table} .= sprintf(
"<g onclick=\"runType(\'TableButton\',evt, \'%s\')\"><path d=\"M%f %f L%f %f A%f %f 0 0 0 %f %f L %f %f A%f %f 0 0 1 %f %f\" stroke=\"black\" fill=\"#bbffbb\" id=\"%s\" data=\"%s\"/>\n",
            $id_in,
            ( $border + $max_radius ) + sin($st) * ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + cos($st) * ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + sin($st) * ( $max_radius - $border ),
            ( $border + $max_radius ) + cos($st) * ( $max_radius - $border ),
            ( $max_radius - $border ),
            ( $max_radius - $border ),
            ( $border + $max_radius ) + sin($en) * ( $max_radius - $border ),
            ( $border + $max_radius ) + cos($en) * ( $max_radius - $border ),
            ( $border + $max_radius ) + sin($en) * ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + cos($en) * ( $max_radius - $border * 1.5 ),
            ( $max_radius - $border * 1.5 ),
            ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + sin($st) * ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + cos($st) * ( $max_radius - $border * 1.5 ),
            "ButtonID" . $table_types[$i],
            $table_types[$i]
        );

        $out{table} .= sprintf(
            "<path d=\"M%f %f A%f %f 0 0 0 %f %f\" fill=\"none\" stroke=\"black\" stroke-width=\"1.5\" id=\"%s\"/>",
            ( $border + $max_radius ) + sin($st) * ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + cos($st) * ( $max_radius - $border * 1.5 ),
            ( $max_radius - $border * 1.5 ),
            ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + sin($en) * ( $max_radius - $border * 1.5 ),
            ( $border + $max_radius ) + cos($en) * ( $max_radius - $border * 1.5 ),
            "ArcButtonID" . $table_types[$i]
        );
        $out{table} .= sprintf(
"<text text-anchor=\"middle\" alignment-baseline=\"middle\"><textPath startOffset=\"%s\" xlink:href=\"#%s\" >%s</textPath></text>\n",
            "50%", "textPath-" . $table_types[$i],
            $table_types[$i]
        );
        $out{table} .= "</g>\n";

    }
    $out{table} .= "</g>\n";

    #Making the actual table-> it starts out as blank
    my $n_r = 7.5;
	$html{table} .= "</svg>\n";
	$html{table} .= sprintf("<div id=\"tableDiv\"  style=\"top:%fpx; position:absolute; left:%fpx; height:%fpx; width:%fpx;\">"
	, $border + $max_radius - 3.5 * $n_r ,  (2.125 * $border),
         $max_radius * 2 - 2.25 * $border,
         $max_radius / 4 + 6 * $n_r );
	   
	$html{table} .= sprintf(
"<table id=\"tableMain\" style=\"position: absolute; left: %fpx; top: %fpx; width: %fpx; height: %fpx; background: grey; font-size:10pt; display: none; table-layout:fixed;\">",
        0,
        4*$n_r,
        $max_radius * 2 - 2.25 * $border,
        $max_radius / 4
    );
    $html{table} .= "</table>";
	$html{table} .= sprintf("<div id=\"searchBox\" style=\"visibility:hidden; position: absolute; left: %fpx; top: %fpx; width: %fpx; height: %fpx; background-color: white; \"><form onsubmit=\"searchInit(event)\" onmouseover=\"showSearch()\" onmouseout=\"searchOff()\" id=\"searchForm\"><input type=\"text\" id=\"searchText\" rows=\"1\" style=\"font-size:14pt;font-color:DarkGray; border: none; height:%s; width:%s;\" onfocus=\"showSearch(\'focus\')\" placeholder=\"Search...\"><input type=\"button\" value=\"Search\" style=\"border: none; background-color: #013220; color: white; padding: 15px 32 px; font-size: 16px;\" onclick=\"searchInit(event)\"></form></div>",
		4*$n_r,
		0,
		 $max_radius * 2 - 2.25 * $border - 10 * $n_r ,
        4* $n_r, 
        4* $n_r,
		$max_radius * 2 - 2.25 * $border - 20 * $n_r . "px"
		
	);
	$html{table} .= sprintf("<svg id=\"searchIcon\" style=\"position: absolute; left: %fpx; top: %fpx; width: %fpx; height: %fpx;\"><rect id=\"searchRect\" x=\"0\" y=\"0\" height=\"%s\" width=\"%s\" style=\"fill:LightGray;\"/ onclick=\"showSearch(\'icon\')\" onmouseover=\"searchOn()\" onmouseout=\"searchOff()\"/><path id=\"magnifyGlass\" transform=\"scale(0.040, 0.04) translate(%s, %s)\" d=\"M754.7 468.7l-22.3-22.3c24.3-33.3 37.5-73.4 37.7-114.7 0-109.5-88.8-198.3-198.3-198.3s-198.3 88.8-198.3 198.3S462.1 530 571.7 530c41.2-.2 81.3-13.4 114.7-37.7l22.3 22.3 152.7 152 45.3-45.3-152-152.6zm-183 0c-75.8 0-137.3-61.5-137.3-137.3S495.8 194 571.7 194 709 255.5 709 331.3c.2 75.7-61 137.1-136.7 137.3-.2.1-.4.1-.6.1z\"></path></svg>", 
		0, 0, 
		3.5* $n_r,
		3.5* $n_r, "100%", "100%", -40 * $n_r, -14 * $n_r
		);
    $html{table} .= sprintf("<a id=\"expandTitle\" xlink:title=\"Expand Table\"><svg id=\"expandSVG\" onclick=\"expandTable()\" style=\"position: absolute; left: %fpx; top: %fpx; width: %fpx; height: %fpx; visibility: hidden;\"><rect id=\"expandRect\" x=\"0\" y=\"0\" height=\"%s\" width=\"%s\" style=\"fill:LightGray;\"/\"/><path id=\"expandArrow\" transform=\"scale(%f, %f)\" d=\"M 25.980469 2.9902344 A 1.0001 1.0001 0 0 0 25.869141 3 L 20 3 A 1.0001 1.0001 0 1 0 20 5 L 23.585938 5 L 13.292969 15.292969 A 1.0001 1.0001 0 1 0 14.707031 16.707031 L 25 6.4140625 L 25 10 A 1.0001 1.0001 0 1 0 27 10 L 27 4.1269531 A 1.0001 1.0001 0 0 0 25.980469 2.9902344 z M 6 7 C 4.9069372 7 4 7.9069372 4 9 L 4 24 C 4 25.093063 4.9069372 26 6 26 L 21 26 C 22.093063 26 23 25.093063 23 24 L 23 14 L 23 11.421875 L 21 13.421875 L 21 16 L 21 24 L 6 24 L 6 9 L 14 9 L 16 9 L 16.578125 9 L 18.578125 7 L 16 7 L 14 7 L 6 7 z\"\"></path></svg></a>", 
		 $max_radius * 2 - 2.25 * $border - 4 * $n_r, 0, 
		3.5* $n_r,
		3.5* $n_r, "100%", "100%", 1.0, 1.0
		);
		
	
    $trans = sprintf(
        "translate(%f,%f) scale(%f, %f)",
       0,
        0,
        0.045, 0.045
    );
    $action = "onclick=\"saveTableTxt()\"";

    #Adding disk image to the table
    $disk_image = $disk_image_orig;
    $disk_image =~ s/TRANS/$trans/g;
    $disk_image =~ s/ACTION/$action/g;
    $disk_image =~ s/TEXT/TXT/g;
    $disk_image =~ s/FS/250/g;
	$saveTableDiskImage = $disk_image;
	
    $html{table} .= sprintf("<div id=\"otherTable\" style=\"position: absolute; left: %fpx; top: %fpx; width: %fpx; height: %fpx; visibility: hidden;\"><svg  width=\"%fpx\" height= \"%fpx\">$saveTableDiskImage</svg></div>", 
	$max_radius * 2 - 2.25 * $border - 4 * $n_r, 0, 0,0, $n_r * 5, $n_r * 4);
	$html{table} .= sprintf("<div id=\"tableType\" style=\"position: absolute; left: %fpx; top: %fpx; width: %fpx; height: %fpx; visibility: hidden; display: inline-block;\" class=\"dropdown\" ><div  class=\"dropbtn\" style=\"width: %fpx; height: %fpx; background-color: LightGray;color: black; text-anchor: middle; alignment-baseline: middle; font-size: 20px; border: none;\" id=\"dropdownID\" onhover=\"showTypeMenu()\">$default_table_type</div><div class=\"drop-content\" id=\"dropContent\">", 
	$max_radius * 2 - 2.25 * $border + 2 * $n_r, 0, 0,0, $n_r * 10, $n_r * 4);
	my $t_cnt = 0;
	foreach my $tableTypes (@table_types)
	{
		if ($tableTypes ne $default_table_type)
		{
			#$html{table} .= sprintf("<button>");
		}
	}
	$html{table} .= "</div></div></div></div>";

    #These are the waiting screens
   	
	$html{end} .= sprintf("<div id=\"Waiting\" style=\"visibility:hidden; top:%fpx; position:absolute;left:%fpx;height:%fpx;width:%fpx;background: rgba(255, 255, 255, 0.50);\" class=\"loader\"></div>",
        ( 2 * $border ),
        1 * $gene_height + ( 2 * $border ),
        $gene_height + ( $border + $max_radius ),
        $gene_height + ( $border + $max_radius )
    );

#$html{end} .= sprintf("<div id=\"WaitingDiv\" style=\"top: %fpx; left: %fpx; position:absolute;height:%fpx;width:%fpx;background: rgba(255, 255, 255, 0.05)\"><h1 id=\"hLoad\" class=\"\" style=\"font-size:500%s; color:red; text-align: center;top: %fpx; left: %fpx; position:absolute;\">Working...</h1></div>",($border + $max_radius)*.6, ($border + $max_radius)*.65,($border + $max_radius)/4, ($border + $max_radius)/4, "%",($border + $max_radius)*.6, ($border + $max_radius)*.65);
    $html{end} .= "</body></html>\n";

    #Getting all the header values associated with the different table types
    my @k = keys(%jso);
    for ( my $i = 0 ; $i < scalar(@k) ; $i++ ) {

        for ( my $j = 0 ; $j < scalar( @{ $jso{ $k[$i] } } ) ; $j++ ) {

            foreach my $kys ( keys( %{ $jso{ $k[$i] }->[$j] } ) ) {

                if ( !$jsHead{ $k[$i] }->{$kys} ) {

                    $jsHead{ $k[$i] }->{$kys} = 1;

                }

            }

        }

    }

    #Writing the table JSON string to the intra-HTML javacript
    #GoodJSON is list of GO and ARO Terms
    $out{script} .= sprintf("\n\nvar goodJSON = \'{");
    for ( my $i = 0 ; $i < scalar(@k) ; $i++ ) {

        $out{script} .= sprintf( "\"%s\":{", $k[$i] );
        for ( my $j = 0 ; $j < scalar( @{ $jso{ $k[$i] } } ) ; $j++ ) {

            if ( $jso{ $k[$i] }->[$j]->{href} ) {

                $out{script} .= sprintf( "\"%s\":1,", $jso{ $k[$i] }->[$j]->{href} );

            }

        }
        chop( $out{script} );
        $out{script} .= "},";

    }
    chop( $out{script} );

#JSON decribing a hash of array of hash giving the table type, the table row, and the header giving the value for the table at that row
    $out{script} .= sprintf("}\';\n\n\nvar tableInfoStr = \'{");
    for ( my $i = 0 ; $i < scalar(@k) ; $i++ ) {

        $out{script} .= sprintf( "\"%s\":[", $k[$i] );
        my @kys = keys( %{ $jsHead{ $k[$i] } } );
        for ( my $j = 0 ; $j < scalar( @{ $jso{ $k[$i] } } ) ; $j++ ) {

            $out{script} .= "{";
            my $cnt1 = 0;
            for ( my $k = 0 ; $k < scalar(@kys) ; $k++ ) {

                if ( $jso{ $k[$i] }->[$j]->{ $kys[$k] } ) {

                    $out{script} .= sprintf( "\"%s\":\"%s\",", $kys[$k], $jso{ $k[$i] }->[$j]->{ $kys[$k] } );
                    $cnt1++;

                }

            }
            if ( $cnt1 > 0 ) {

                chop( $out{script} );

            }
            $out{script} .= "},";

        }
        if ( scalar(@kys) > 0 ) {

            chop( $out{script} );

        }
        $out{script} .= "],";

    }
    chop( $out{script} );
    $out{script} .= "}\';\n";

    #JSON string giving the hash of the array giving the column name for each of the columns for each table type
    $out{script} .= "var tableHeadStr = \'{";
    for ( my $i = 0 ; $i < scalar(@k) ; $i++ ) {

        $out{script} .= sprintf( "\"%s\":[", $k[$i] );
        for ( my $k = 0 ; $k < scalar( @{ $head{ $k[$i] } } ) ; $k++ ) {

            $out{script} .= sprintf( "\"%s\",", $head{ $k[$i] }->[$k] );

        }
        chop( $out{script} );
        $out{script} .= "],";

    }
    chop( $out{script} );

    $out{script} .= "}\';\n ";

    #JSON string giving the hash of the array giving the column wdith for each of the columns for each table type

    $out{script} .= "var tableWidthStr = \'{";
    for ( my $i = 0 ; $i < scalar(@k) ; $i++ ) {

        $out{script} .= sprintf( "\"%s\":[", $k[$i] );
        for ( my $k = 0 ; $k < scalar( @{ $col_width{ $k[$i] } } ) ; $k++ ) {

            $out{script} .= sprintf( "\"%s\",", ( $col_width{ $k[$i] }->[$k] * 100 )  );

        }
        chop( $out{script} );
        $out{script} .= "],";

    }
    chop( $out{script} );

    #JSON giving a hash that mapes each gene to a region
    $out{script} .= "}\';\n ";
    $out{script} .= "var region_2_gene = \'{";
    my @k2 = keys(%reg2gene);
    for ( my $i = 0 ; $i < scalar(@k2) ; $i++ ) {

        $out{script} .= sprintf( "\"%s\":\"%s\",", $k2[$i], $reg2gene{ $k2[$i] } );

    }
    chop( $out{script} );
    $out{script} .= "}\';\n ";
    $out{script} .= "var term_2_type = \'{";
    @k2 = keys(%term2type);
    for ( my $i = 0 ; $i < scalar(@k2) ; $i++ ) {

        $out{script} .= sprintf( "\"%s\":\"%s\",", $k2[$i], $term2type{ $k2[$i] } );

    }
    chop( $out{script} );
    $out{script} .= "}\';\n ";

    $out{script} .= "var type_2_term = \'{";
    @k2 = keys(%type2term);
    for ( my $i = 0 ; $i < scalar(@k2) ; $i++ ) {

        $out{script} .= sprintf( "\"%s\":[", $k2[$i] );
        for ( my $j = 0 ; $j < scalar( @{ $type2term{ $k2[$i] } } ) ; $j++ ) {

            $out{script} .= sprintf( "\"%s\",", $type2term{ $k2[$i] }->[$j] );

        }
        chop( $out{script} );
        $out{script} .= "],";

    }
    chop( $out{script} );
    $out{script} .= "}\';\n ";

    #intra-html javascript is finished and the main.js script is also added
    $out{script} .= "
	//]]>
	</script><script type=\"text/javascript\" src=\"scripts/main.javascript.js\"></script>\n";

    #end of the beginning string of the main page HTML
    $html{beg} .=
      "</head><body xmlns=\"http://www.w3.org/1999/xhtml\"  xmlns:xlink=\"http://www.w3.org/1999/xlink\" onload=\"init()\">";

    #writing HTML opening lines, loading screen, intra-main page javascript, and the end of the html page
    print SVG_MAIN $html{beg};
    print SVG_MAIN $out{load};
    print SVG_MAIN $out{script};
    print SVG_MAIN $out{beg};

    #Making all the thumbnail sketches for the left side menu
    my $thumb_cnt   = 0;     #Number of thumbnails in the menu
    my $type_cnt    = 0;     #number of types in the menu
    my $type_height = 25;    #
    my $max_chr     = 0;

    #removed white-space: nobreak; justify-content: top;
    #This is the menu bar that is shown when the menu is out i.e. visible
    $out{thumb}->{beg} .= sprintf(
"<div id = \"coreMenuOut\" style=\"background-color: #bbbbbb; top: %fpx; width: %fpx; height: %fpx; left: 0px; position: absolute; padding: 0; margin: 0; font-size: 20px; transform-origin: 0 100%s; transform: rotate(90deg);  text-align: center;\" onmouseover=\"turnOnHighlight(event)\" onmouseout=\"turnOffHighlight(event)\" onclick=\"selectMenu(event)\"> Assembly Core Regions</div>",
        -1 * $gene_height,
        ( 2 * ( $border + $max_radius ) ),
        1 * $gene_height,
        "%",
        ( 2 * ( $border + $max_radius ) ),
        1 * $gene_height
    );

    #This is the menu bar that is shown when the menu is in i.e. not visible
    $out{thumb}->{beg} .= sprintf(
"<div id = \"coreMenuIn\" style=\"background-color: #dddddd; top:0px; height: %fpx; width: 0px; left: 0px; position: absolute; overflow-y:auto;\">",
        ( 2 * ( $border + $max_radius ) ),
        1 * $gene_height,
        ( 2 * ( $border + $max_radius ) ),
        1 * $gene_height
    );
    my $chr_type_list;

    #Cycling through the chromosome types
    for my $thumb_type ( keys(%thumb_list) ) {

        #For each chromosome type, make a menu header that can be expanded or hidden
        $out{thumb}->{$thumb_type} .= sprintf(
"<p Showing=\"on\" id =\"$thumb_type\" type=\"$thumb_type\" onclick=\"changeMenu(event)\" onmouseover=\"turnOnHighlight(event)\" onmouseout=\"turnOffHighlight(event)\" style=\"text-align:center; font-size: 20px; background-color: #bbbbbb;\">%s</p>\n",
            $thumb_type . "(-)" );

        #Going through all the Thumbnails of that chromosome type and appends it to the svg
        for my $cur_chr ( @{ $thumb_list{$thumb_type} } ) {

            #saving the thumbnail type in the chromosome type array
            $chr_type_list->[$cur_chr] = $thumb_type;
            if ( $cur_chr > $max_chr ) { $max_chr = $cur_chr; }
            my $id_in = "svgOut" . $cur_chr;
            my $vis   = "hidden";              #All chromsomes are oringally hidden; on the initialization, the lowest

            #Writing the thumbnail SVG and the main SVG
            print SVG_MAIN "<svg id=\"$id_in\" visibility=\"$vis\">";

            $out{thumb}->{$cur_chr} .= sprintf(
"<svg id=\"svgThumb$cur_chr\" type=\"$thumb_type\" visibility=\"hidden\" viewBox=\"0 0 %f %f\" x=\"0\" y=\"%f\" width=\"%f\" height=\"%f\" onclick=\"setSVG($cur_chr)\">",
                $thumb_size_w, $thumb_size_h, ($thumb_cnt) * $thumb_size_h,
                $thumb_size_w, $thumb_size_h
            );

            #shrinking the thumbnail to predefined sized
            $out{thumb}->{$cur_chr} .= "<g transform=\"scale($thumb_dif,$thumb_dif)\">";

            #prints the chromosome number in the middle of the thumbnail
            $out{thumb}->{$cur_chr} .= sprintf(
                "<text x=\"%f\" y=\"%f\" font-size=\"%f\" alignment-baseline=\"middle\" text-anchor=\"middle\">%s</text>",
                ( $border + $max_radius ),
                ( $border + $max_radius ),
                ( $border + $max_radius ) * 2 / 3, $cur_chr
            );
            for ( my $i = 100 ; $i >= 2 ; $i-- ) {

                if ( $out{$cur_chr}->{$i} ) {

                    print SVG_MAIN $out{$cur_chr}->{$i};
                    if ( $i == 100 ) {

                        #thumbnail only has level 1 and level 100
                        #also all the mouseevents are turned off
                        my $new = $out{$cur_chr}->{$i};
                        $new =~ s/onclick/oldonclick/g;
                        $new =~ s/onmouseout/oldmouseout/g;
                        $new =~ s/onmouseover/oldonmouseover/g;

                        #writing these levels to the thumbnail SVG
                        $out{thumb}->{$cur_chr} .= $new;

                    }

                }

            }

            #closing the thumbnail SVG
            print SVG_MAIN "</svg>";
            $thumb_cnt++;
            $id_in = "svgPreview" . $cur_chr;

            #writing all the preview panes of the main chromosome to the SVG
            print SVG_MAIN "<svg id=\"$id_in\" visibility=\"$vis\">";
            for ( my $i = 100 ; $i >= 0 ; $i-- ) {

                if ( $out{$cur_chr}->{preview}->{$i} ) {

                    print SVG_MAIN $out{$cur_chr}->{preview}->{$i};

                }

            }

            #Closing the svg
            print SVG_MAIN "</svg>";

            #writing the levels (0-1) of the main chromosome to the SVG
            $id_in = "svgArc" . $cur_chr;
            print SVG_MAIN "<svg id=\"$id_in\" visibility=\"$vis\">";

            for ( my $i = 1 ; $i >= 0 ; $i-- ) {

                if ( $out{$cur_chr}->{$i} ) {

                    print SVG_MAIN $out{$cur_chr}->{$i};
                    if ( $i == 1 || $i == 100 ) {

                        my $new = $out{$cur_chr}->{$i};
                        $new =~ s/onclick/oldonclick/g;
                        $new =~ s/onmouseout/oldmouseout/g;
                        $new =~ s/onmouseover/oldonmouseover/g;
                        $out{thumb}->{$cur_chr} .= $new;

                    }

                }

            }
            $out{thumb}->{$cur_chr} .= "</g></svg>";

            print SVG_MAIN "</svg>";

        }

    }
    if ( ($thumb_cnt) * $thumb_size_h < ( 2 * ( $border + $max_radius ) ) ) {

        $out{thumb}->{end} .= sprintf(
"<svg x=\"0\" y=\"%f\" id=\"blankSVG\" width=\"%f\" height=\"%f\" visibility=\"hidden\" background=\"red\" onclick=\"runType(\'selectMenu\',event, \'\')\"></svg>",
            ( ($thumb_cnt) * $thumb_size_h ),
            $thumb_size_w, ( 2 * ( $border + $max_radius ) ) - ( 25 + ($thumb_cnt) * $thumb_size_h )
        );

    }
    $out{thumb}->{end} .= "</div>";

    #Making the Legend

    print SVG_MAIN $out{table};
    print SVG_MAIN $html{table};

    print SVG_MAIN $out{thumb}->{beg};

    #Writing all the thumbnails on the SVG
    for my $thumb_type ( "Assembly_core", "Cycle fGRs", "Other fGRs" ) {

        if ( $out{thumb}->{$thumb_type} ) {

            print SVG_MAIN $out{thumb}->{$thumb_type};
            for ( my $cur_chr = 1 ; $cur_chr <= $max_chr ; $cur_chr++ ) {

                if ( $chr_type_list->[$cur_chr] eq $thumb_type ) {

                    print SVG_MAIN $out{thumb}->{$cur_chr};

                }

            }

        }

    }
    print SVG_MAIN $out{thumb}->{end};

    print SVG_MAIN $out{end};
    print SVG_MAIN $html{end};

}

#Creating a new string that contains the entire FGR page (without javascript which is seperately written)
#The string is an SVG
sub make_fgi_page_svg {

    my ( $file, $data, $ont, $out ) = @_;
    my $height  = scalar( @{ $data->{order} } ) + 11.5;
    my $sz      = $gene_height * 1.0;                     #Setting the size of the genes in the arrow to the gene height
    my $sz_h    = $sz;                                    #height size
    my $tot     = 0;                                      #This is the total width size with equal sized genes
    my $tot_old = 0;                                      #This is the old total width size. Not used
    my %legend;                                           #hash that matches function name with the color in the legend
    my $sz_dif     = 0.05;                                #
    my $max_st     = 75 / $sz_dif;                        #This is the default size of the boundary regions
    my $max_end    = 75 / $sz_dif;                        #This is the default size of the boundary regions
    my $std_gene_x = $max_st / 2;                         #Standard gene size- used in the size
    my $gene_list  = "";                                  #keep a list of the genes in the

    ###################################################################################
    #Step (1): Going through the genes to (a) get the appropiate width of the svg and (b) set the colors + ids of the legend
    ###################################################################################

    for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

        #Getting the gene information for every gene in the region
        my ( $ID, $num, $type, $f_id, $member, $cl_id, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[$j] };
        $tot_old += abs( $end_loc - $st_loc );    #The old version where the size is based on the gene size
        $tot += $std_gene_x;                      #Getting the total horizontal size of the
        my $new_c = $c;                           # color of the gene
        $new_c =~ s/\s//g;                        #removing any spaces from the color
        if ( !$new_c ) { $new_c = "000000"; }     #setting default value to black
        $legend{ $info->{type_ref} } = "#" . $new_c;    #setting the color of on the legend
        $gene_list .= "\"$cl_id\":\"#$new_c\",";        #add it

    }
    chop($gene_list);                                   #remove the comma from the end

    my $tot_width = max( $tot * $sz_dif, 120 );         #getting the size of the

    #Going through the boundary regions-> resetting the maximium boundary gene sizes as needed
    #Also adding the function values to the legend as needed
    foreach my $j ( keys( %{ $data->{st_sz} } ) ) {

        if ( $max_st < $data->{st_sz}->{$j} ) {

            #$max_st = $data->{st_sz}->{$j};
        }
        if ( $max_end < $data->{end_sz}->{$j} ) {

            #$max_end = $data->{end_sz}->{$j};
        }
        if ( $data->{end_type}->{$j} ) {

            $legend{ $data->{end_type}->{$j} } = "#" . $data->{end_col}->{$j};

        }
        if ( $data->{st_type}->{$j} ) {

            $legend{ $data->{st_type}->{$j} } = "#" . $data->{st_col}->{$j};

        }

    }

    #Adding the bounding genes sizes to the image width size (both old and new)
    $tot_old += $max_end + $max_st;
    $tot     += $max_end + $max_st;

    ###################################################################################
    #Step (2) Writing the gene and genome information in JSON
    ###################################################################################

    #Starting the intrahtml javascript
    my $jscript = "<script type=\"text/javascript\">
	//<![CDATA[

	var currentID=\"\";";
    my $terms;

    $gene_list =~ s/\n//g;

    #adding the list of genes and their colors to the javascript
    $jscript .= sprintf("\n\nvar geneListJSON = \'{$gene_list}\';");

    #Adding the names in order to a javascript array JSON string
    $jscript .= sprintf("\n\nvar namesJSON = \'[");
    for ( my $j = 0 ; $j < scalar( @{ $data->{order} } ) ; $j++ ) {

        my $name = $data->{order}->[$j];
        $jscript .= sprintf( "\"%s\",", $name );

    }
    chop($jscript);
    $jscript .= "]\';\n";

    $jscript .= "var isGenomeOn = new Array();\n";
    for ( my $j = 0 ; $j < scalar( @{ $data->{order} } ) ; $j++ ) {

        my $name    = $data->{order}->[$j];
        my $genomes = $data->{gen}->{$name};
        my @genList = split ":", $genomes;
        foreach my $a (@genList) {

            $jscript .= sprintf( "isGenomeOn[\'%s\']=0; ", $a );

        }

    }
    $jscript .= "\n";

#Writes the JSON string that contains the starting bounding gene and the ending bounding gene and counts for each of the genes
    my $startsJSON = "var startJSON=\'{";
    my $endsJSON   = "var endJSON=\'{";
    my $cntsJSON   = "var cntJSON=\'{";
    $jscript .= sprintf("\n\nvar genomeJSON = \'{");
    for ( my $j = 0 ; $j < scalar( @{ $data->{order} } ) ; $j++ ) {

        my $name    = $data->{order}->[$j];
        my $genomes = $data->{gen}->{$name};
        my @genList = split ":", $genomes;
        $startsJSON .= sprintf( "\"%s\":\"%s\",", $name, $data->{st_id}->{$name} );
        $endsJSON   .= sprintf( "\"%s\":\"%s\",", $name, $data->{end_id}->{$name} );
        $cntsJSON   .= sprintf( "\"%s\":\"%s\",", $name, $data->{cnt}->{$name} );

        $jscript .= sprintf( "\"%s\":[", $name );
        foreach my $a (@genList) {

            $jscript .= sprintf( "\"%s\",", $a );

        }
        chop($jscript);
        $jscript .= "],";

    }
    chop($jscript);
    $jscript .= "}\';\n";
    chop($startsJSON);
    chop($endsJSON);
    chop($cntsJSON);

    #adding all to the javascript string
    $jscript .= $startsJSON . "}\';\n" . $endsJSON . "}\';\n" . $cntsJSON . "}\';\n";

    #Making a 2nd JSON string for which every gene points to the genomes that contain this gene
    $jscript .= sprintf("\n\nvar geneJSON = \'{");
    for ( my $j = 0 ; $j < scalar( @{ $data->{order} } ) ; $j++ ) {

        my $name    = $data->{order}->[$j];
        my $genomes = $data->{gen}->{$name};
        my @genList = split ":", $genomes;
        $jscript .= sprintf( "\"%s\":[", $name );
        foreach my $a (@genList) {

            $jscript .= sprintf( "\"%s\",", $a );

        }
        chop($jscript);
        $jscript .= "],";

    }
    chop($jscript);
    $jscript .= "}\';\n";

    #Making the JSON string necessary to make the gene main information page
    $jscript .= sprintf("\n\nvar fastaJSON = \'{");
    for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

        my ( $ID, $num, $type, $f_id, $member, $cl_id, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[$j] };
        $cl_id =~ /_(\d+)/;
        my $a             = $1;
        my $genomeList    = "";
        my $genomeLineLen = 0;

        foreach my $genomeID ( keys( %{ $seqs{byClust}->{$a} } ) ) {

            $genomeLineLen += length($genomeID);
            if ( $genomeLineLen > 80 ) { $genomeList .= "<br>"; $genomeLineLen = length($genomeID); }
            $genomeList .= "$genomeID,";

        }
        if ($genomeList) {

            chop($genomeList);

        }

        $jscript .= sprintf(
"\"%s\":{\"Region\":\"%s\",\"Cluster\":\"%s\", \"Name\":\"%s\",\"# of Genomes\":\"%s\",\"Functional IDs\":\"%s\",\"Associated Terms\":\"%s\",\"Sequence\":\"%s\", \"Mean\":\"%s\",\"Standard Deviation\":\"%s\",\"Maximum Length\":\"%s\",\"Minimum Length\":\"%s\", \"Genomes\":\"%s\"},",
            $cl_id, "FGI_" . $f_id,
            $cl_id, $type, $num, $info->{type_ref}, $info->{oth_ref}, $seqs{cent}->{$a}->{seq},
            $info->{Mean},
            $info->{"Standard Deviation"},
            $info->{"Maximum Length"},
            $info->{"Minimum Length"}, $genomeList
        );
        if ( $info->{oth_ref} ) {

            my @term_list = split ";", $info->{oth_ref};
            foreach my $a (@term_list) {

                if ( $ont{$a} ) { $terms->{$a} = $ont{$a}; }

            }

        }

    }
    chop($jscript);
    $jscript .= "}\';\n";

  #Writing the sequence/fasta for every gene in a cluster information to a new JSON file: this speeds up the loading process.
    my $tmp_fasta;    #This stores the
    for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

        #ID = gene, $num = , type = gene_type, f_id = , member
        my ( $ID, $num, $type, $f_id, $member, $cl_id, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[$j] };
        $cl_id =~ /_(\d+)/;
        my $a = $1;
        my $nJscript .= sprintf( "var currentID =\'%s\';\n\nvar allFastaJSON = \'[", $cl_id );
        my $nJscriptAdd = "";
        foreach my $genomeID ( keys( %{ $seqs{byClust}->{$a} } ) ) {

            $nJscriptAdd .= sprintf(
                "{\"seq\":\"%s\", \"id\":\"%s\"},",
                $seqs{byClust}->{$a}->{$genomeID},
                $seqs{byClustId}->{$a}->{$genomeID}
            );
            if ( $st_loc < $end_loc ) { $tmp_fasta->{$genomeID} .= $seqs{byClust}->{$a}->{$genomeID}; }
            else {

                my $tmp_seq = reverse( $seqs{byClust}->{$a}->{$genomeID} );
                if ($tmp_seq) { $tmp_seq =~ tr/ATCGatcg/TAGCtagc/; }
                $tmp_fasta->{$genomeID} .= $tmp_seq;

            }

        }
        if ($nJscriptAdd) {

            chop($nJscriptAdd);

        }
        $nJscript .= $nJscriptAdd . "]\';\n";

        #Writing the gene cluster sequence JSON to the file allFasta.json
        open( OUTJSON, ">", $out_dir . "/$output_id/json/" . ($cl_id) . ".allFasta.json" );
        print OUTJSON $nJscript;
        close(OUTJSON);

    }

    #Writing the cluster centroid sequence to a JSON string
    $jscript .= sprintf("\n\nvar seqsJSON = \'{");
    foreach my $id1 ( keys(%$tmp_fasta) ) {

        $jscript .= sprintf( "\"%s\":\"%s\",", $id1, $tmp_fasta->{$id1} );

    }
    chop($jscript);
    $jscript .= sprintf("}\';\n");

    #Writing the cluster GO Terms to a JSON string
    $jscript .= sprintf("\n\nvar termJSON = \'{");
    my $add_on = "";    #used to make sure that there is something in the JSON string
    foreach my $a ( keys(%$terms) ) { $add_on .= sprintf( "\"%s\":\"%s\",", $a, $terms->{$a} ); }
    if ($add_on) { chop($add_on); }
    $jscript .= $add_on . "}\';\n";

    $jscript .= "//]]>";

    #Add the tree JSON string as needed
    if ( $tree_level > 0 ) {

        $jscript .= "</script><script type=\"text/javascript\" src=\"../json/tree.json\">";

    }

    #Finish the javascript variable assignment (all JSON strings)

    $jscript .= "
</script><script type=\"text/javascript\" src=\"../scripts/fgi.functions.js\"></script>\n";

    ###################################################################################
    #Step (3) Writing the legend image
    ###################################################################################

    my @leg_id = keys(%legend);    #All the functions in the legend
    my $x1     = 0;                #location on the x-axis
    my $y_len  = 1;                #hieght of this particular rows
    my @max_n;                     #array that stores the Maximum number of text lines in each of the horizontal rows
    my @x_len;                     #array that stores the x location of each legend key
    my @y_loc;                     #array that stores the y location of each legend key
    my $x_l = 0;                   #legend text horizontal size. Reset at each number
    my @y_l;                       #array that keeps track of the y location os
    my $cur_y         = 0;              #The current vertical level of the legend
    my $legend_width  = 0;              # the width of the legend on the bottom
    my $legend_height = 0;              #The height of the legend
    my $sz_h2         = $sz_h * 0.5;    #Another gene height that is set at half the height

    #Getting the x and y location on the legend for each of the functions in the fGI
    #Also helping to set up the svg
    my $r1 = 0.25;    #Horizontal size of a letter
    $max_n[ $y_len - 1 ] = 0;    #
    for ( my $i = 0 ; $i < scalar(@leg_id) ; $i++ ) {

        my $y1      = $sz_h * ( $height + 1.5 ) + 2;                     #The vertical location
        my $max_l   = 0;                                                 #Maximum vertival size in this row
        my $id_name = split_id( $leg_id[$i], $x1 + $sz_h * 0.6, 15 );    #divides an id into different lines as needed
        my $n       = -1;
        while ( $id_name =~ /([^\>\<]+)\<\/tspan\>/g )                   #Getting the number of lines in a in id
        {

            $n++;
            if ( length($1) > $max_l ) { $max_l = length($1); }          # re

        }

        $x_l = $max_l * $sz_h * $r1 + $sz_h * 0.75;
        if ( $x1 + $x_l >
            $tot_width )  #if the legend exceeds the horizontal size, this moves down to lower y-value and resets the x-value
        {

            $x_len[ $y_len - 1 ] = $x1;    #resetting x-value to the start
            $legend_height += ( 1 + $max_n[ $y_len - 1 ] );            #adding to the legend height
            $cur_y += $sz_h * ( ( $max_n[ $y_len - 1 ] + 1 ) * 1 );    #
            if ( $x1 + $x_l > $legend_width ) { $legend_width = $x1 + $x_l + 4; }    #resetting the legend width as needed
            $y_len++;
            $x1 = 0;

        }
        if (  !$max_n[ $y_len - 1 ]
            || $n >
            $max_n[ $y_len - 1 ] )    #if the number in a line exceeds the current maximium, makes thje max equal this number
        {

            $max_n[ $y_len - 1 ] = $n;

        }
        $x1 += $x_l;                  #Adding to x value
        $y_loc[$i] = $y_len - 1;      #Keeping y as needed
        $y_l[$i]   = $cur_y;

    }

    $x_len[ $y_len - 1 ] = $x1;
    $legend_height += $cur_y + ( $max_n[ $y_len - 1 ] + 1 );
    if ( $x1 > $legend_width ) { $legend_width = $x1 + 4; }

    #Starting the FGI html page
    $out = sprintf(
"<html><body xmlns=\"http://www.w3.org/1999/xhtml\"  xmlns:xlink=\"http://www.w3.org/1999/xlink\" onload=\"makeImage()\" onresize=\"makeImage()\">"
    );

    #Adding the javascript
    $out .= $jscript;

    #Starting to write the legend using the locations generated above

    #Setting the x location of the legend box
    $x1 = ( $tot_width * $sz_dif - $legend_width ) / 2;
    if ( $x1 < 0 ) {

        $x1           = 0;
        $legend_width = $tot_width * $sz_dif;

    }    #If x1 is greater than the total width, then just set it to zero
    $legend_width = max( 250, $legend_width );
    $out .= sprintf(
"<div id = \"legendDiv\" ht=\"%f\" wd=\"%f\"  position=\"fixed\" num=\"%f\"><svg id=\"legendSVG\" version=\"1.2\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id=\"$file\" height=\"\%f\" width=\"\%f\">\n",
        ( $legend_height + 1.5 ) * $sz_h,
        $legend_width, scalar(@leg_id), ( $legend_height + 1.5 ) * $sz_h,
        $legend_width
    );
    $out .= sprintf("<g id=\"legend\">\n");

    #Adding background rectangle
    $out .=
      sprintf( "<rect id=\"legRect\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"#eeeeee\" stroke =\"black\"/>",
        $x1, 0, $legend_width, ( $legend_height + 1.5 ) * $sz_h );

    #Adding the text title
    $out .= sprintf(
        "<text id=\"LegText\" x=\"%f\" y=\"%f\" text-anchor=\"middle\" font-size=\"%s\">Legend</text>\n",
        $x1 + $legend_width / 2,
        $sz_h * ( $legend_height + 1.25 ) + 2,
        $sz_h * 1.25
    );

    $y_len = 0;

    #Writing all the functional annotation
    for ( my $i = 0 ; $i < scalar(@leg_id) ; $i++ ) {

        my $x_l     = 0;
        my $id_name = split_id( $leg_id[$i], $x1 + $sz_h * 0.6, 15 );    #the functional name as a multiple line
        my $max_l   = 0;
        my $n       = -1;
        while ( $id_name =~ /([^\>\<]+)\<\/tspan\>/g ) {

            $n++;
            if ( length($1) > $max_l ) { $max_l = length($1); }

        }

        $x_l = $max_l * $sz_h * $r1 + $sz_h * 0.75;                      #getting the width of the legend text

        if ( $y_loc[$i] != $y_len ) { $y_len++; $x1 = ( $tot_width - $x_len[$y_len] ) / 2; }
        my $y1 = 0;
        if ( $y_l[$i] ) { $y1 = $y_l[$i]; }
        $id_name = split_id( $leg_id[$i], $sz_h * 0.6, 15 );             #resetting the id name by removing the
               #Writing out the legend functional text and colored square in an SVG
        $out .= sprintf( "<svg id=\"%s\"  x=\"%f\" y=\"%f\">", "TextID" . $i, $x1, $y1 );
        $out .= sprintf(
            "<rect  href=\"%s\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" />\n",
            $i, 0,
            $sz_h * ( $max_n[$y_len] ) / 5,
            $sz_h * 0.5,
            $sz_h / 3.0,
            $legend{ $leg_id[$i] }
        );
        $out .= sprintf(
"<text href=\"%s\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"auto\" font-size=\"%f\" dominant-baseline=\"hanging\" text-anchor=\"start\" fill=\"black\">%s</text>\n",
            $i,
            $sz_h * 0.5,
            $sz_h * ( $max_n[$y_len] - $n ) / 5,
            $border * 1.15,
            $sz_h * .5, $id_name
        );
        $out .= "</svg>";
        $x1 += $x_l;

    }
    $out .= "</g></div>\n";

    ###################################################################################
    #Step (4) Making the functional buttons at the bottom of the screen that can allows
    ###################################################################################

    #Drawing the buttons on the left side of the legend. These are the save
    #Drawing a disk image for the saving the image
    my $trans      = sprintf( "translate(%f,%f) scale(%f, %f)", 0, 10, 0.045, 0.045 );
    my $action     = "onclick=\"saveFullSVG(0,\'$file\')\"";
    my $disk_image = $disk_image_orig;
    $disk_image =~ s/TRANS/$trans/g;
    $disk_image =~ s/ACTION/$action/g;
    $disk_image =~ s/TEXT/FULL/g;
    $disk_image =~ s/FS/250/g;

    $out .= sprintf( "<div id = \"saveDiv\" position=\"absolute\" ht=\"%f\" wd=\"%f\">", $sz_h * 2 + 120, $sz_h * 3 + 10 );
    $out .= sprintf(
"<svg id=\"diskSVG\" version=\"1.2\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
    );
    $out .= $disk_image;    #This is save the window image disk

    $trans      = sprintf( "translate(%f,%f) scale(%f, %f)", 0, 45, 0.045, 0.045 );
    $action     = "onclick=\"saveFullSVG(1,\'$file\' )\"";
    $disk_image = $disk_image_orig;
    $disk_image =~ s/TRANS/$trans/g;
    $disk_image =~ s/ACTION/$action/g;
    $disk_image =~ s/TEXT/WINDOW/g;
    $disk_image =~ s/FS/150/g;

    $out .= $disk_image;    #This is save the full image disk
                            #Drawing the "slider" to change whether the image is saved as a PNG or SVG
    $out .= sprintf(
"<g id=\"ChangeSVG\" onclick=\"changeOut(\'svg\')\"><rect id=\"rectsvg\" y=\"%f\" x=\"0\" width=\"%f\" height=\"%f\" fill-opacity=\"0.5\" fill=\"#bb00bb\"/><text y=\"%f\" x=\"%f\">SVG</text></g>",
        85, $sz_h * 3, $sz_h, 82 + $sz_h, 10 );
    $out .= sprintf(
"<g id=\"ChangePNG\" onclick=\"changeOut(\'png\')\"><rect id=\"rectpng\" y=\"%f\" x=\"0\" width=\"%f\" height=\"%f\" fill-opacity=\"0.5\" fill=\"#bbbbbb\"/><text y=\"%f\" x=\"%f\">PNG</text></g>",
        85 + $sz_h, $sz_h * 3, $sz_h, 82 + 2 * $sz_h, 10 );

    $out .= "</svg></div>";

#The right bottom set of buttons is in a 2x2 set: trimming the fGR, remove/show the singleton, saving the sequences as a fasta of the save fGIs,
#All of these call the javascript
    $out .= sprintf("<div id=\"sortDiv\">");
    $out .= sprintf(
"<button style=\"left:1; top:5; position: absolute; width: %f; height: %f; background-color:#4444ff; \" onclick=\"trimGenomesImage()\" trim=\"0\" id=\"trimButton\">Trim Rows</button>",
        $sz_h * 3, $sz_h * 2, $sz_h / 4, $sz_h / 4 );
    $out .= sprintf(
"<button style=\"left:1; top:65;position: absolute; width: %f; height: %f; background-color:#bb00bb; \" onclick=\"showSingleton()\" trim=\"0\" id=\"singletonButton\">Remove Singletons</button>",
        $sz_h * 3, $sz_h * 2, $sz_h / 4, $sz_h / 4 );
    $out .= sprintf(
"<button style=\"left:%f; top:5;position: absolute; width: %f; height: %f; background-color:#00bb00\; \" onclick=\"saveMultiFasta()\" id=\"saveFastaButton\">Save Fasta</button>",
        $sz_h * 3.25,
        $sz_h * 3, $sz_h * 2, $sz_h / 4, $sz_h / 4
    );

    #Only has the draw phylogeny button if a tree is loaded
    if ( $tree_level > 0 ) {

        $out .= sprintf(
"<button style=\"left:%f; top:65;position: absolute; width: %f; height: %f; background-color:#00bbbb\; \" onclick=\"loadTreePage()\" id=\"loadTreeButton\">View Tree</button>",
            $sz_h * 3.25,
            $sz_h * 3, $sz_h * 2, $sz_h / 4, $sz_h / 4
        );

    }
    $out .= "</div>";

    $tot -= $max_end - $max_st;

    ###################################################################################
    #Step (5) Writing the bounding images, functions and button to out
    ###################################################################################

    my $cnt     = 0;     #x location in the svg
    my $cnt_old = 0;
    my $cur_fam = "";    #The current fGI: also used for the background color: the color is changed at a switch of fGI
    my $i  = scalar( @{ $data->{order} } );    # The number of genes in the fGI (not counting the bounding genes)
    my $i2 = 0;
    my $even = 0;    #Whether the column is even or odd in order to do the zebra background for the fGI
    my @back_col    = ( "#eeeeee", "#cccccc" );    #this is the alternating vertical background zebra stripe colors
    my $fam_st      = $cnt;                        #The starting x-value of the fGI: used to draw the background
    my $fam_st_old  = $cnt_old;                    #This is the same but the old version: ie when genes were sized by # of bp
    my $fam_end     = 0;                           #The ending x-value of the fGI: used to draw the background
    my $fam_end_old = 0;                           #This is the same but the old version: ie when genes were sized by # of bp
    my $dif_i       = 4;
    my $last_num    = 0;

    #total size of the fGR window
    my $tot_h = $sz_h * ( $height + $dif_i ) + $legend_height;
    my $tot_w = max( 260, 10 + $tot * $sz_dif );

    my $after = sprintf(
"<div id = \"mainDiv\" tot_width=\"$tot_w\" tot_height=\"$tot_h\"  position=\"fixed\" onscroll=\"changeWithScroll()\">"
    );
    $after .= sprintf(
"<svg id=\"mainSVG\" version=\"1.2\"  baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" style=\"top:%f; left:0; height:%f; width:%f; position: absolute;\"id=\"$file\" numRow=\"%d\">\n",
        $sz_h * $dif_i,
        $sz_h * ( scalar( @{ $data->{order} } ) ),
        ($tot_width), scalar( @{ $data->{order} } )
    );

    #This is the top row: ie the name of the fGI
    my $top_row = sprintf( "<svg id=\"topRow\" style=\"top:0; left:0; height:%f; width:%f; position: absolute;\">",
        $sz_h * $dif_i, $tot_width );

    #This is the top row: ie the name of the gene cluster
    my $bottom_row = sprintf(
        "<svg id=\"bottomRow\" style=\"top:%f; left:0; height:%f; width:%f; position: absolute;\">",
        $sz_h * ( $dif_i + scalar( @{ $data->{order} } ) ),
        ( 4.75 * $sz_h ), $tot_width
    );

    #Titles of the two
    $bottom_row .= sprintf(
        "<text id=\"footerText\" x=\"%f\" y=\"%f\" text-anchor=\"middle\" font-size=\"%s\">Cluster ID</text>\n",
        ( $sz_dif * ($tot) + 10 ) / 2,
        $sz_h * (4.5), $sz_h
    );
    $top_row .= sprintf(
        "<text id=\"headerText\" x=\"%f\" y=\"%f\" text-anchor=\"middle\" font-size=\"%s\">fGI ID</text>\n",
        ( $sz_dif * ($tot) + 10 ) / 2,
        $sz_h * 1.35, $sz_h
    );

    #Going through the genes in order to draw the background
    for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

        my ( $ID, $num, $type, $f_id, $member, $geneID, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[$j] };
        my $gene_size_old = abs( $end_loc - $st_loc );  #Old version: gene size is the bp
        my $gene_size     = $std_gene_x;                #This is the gene sized used- a default size where each gene is equal

        $fam_st     = $cnt;                             #the current location (new style)
        $fam_st_old = $cnt_old;                         #current location (old style)
        if ( $cur_fam && $cur_fam ne $f_id ) {

            $even = ( $even == 0 );                     #setting the background color

        }

        #drawing the backround color in the main, top and bottom
        $after .= sprintf(
"<rect id=\"%s\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" old_x=\"%f\" old_width=\"%f\" type=\"background\"/>\n",
            "main" . $geneID,
            ($fam_st) * $sz_dif,
            (0) * ($sz_h),
            ($gene_size) * $sz_dif,
            ( $i - $i2 ) * ($sz_h),
            $back_col[$even],
            ($fam_st_old) * $sz_dif,
            ( $fam_end_old - $fam_st_old ) * $sz_dif
        );
        $top_row .= sprintf(
            "<rect id=\"%s\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" old_x=\"%f\" old_width=\"%f\"/>\n",
            "top" . $geneID,
            ($fam_st) * $sz_dif,
            (1.75) * ($sz_h),
            ($gene_size) * $sz_dif,
            ($dif_i) * ($sz_h),
            $back_col[$even],
            ($fam_st_old) * $sz_dif,
            ( $fam_end_old - $fam_st_old ) * $sz_dif
        );
        $bottom_row .= sprintf(
            "<rect id=\"%s\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" old_x=\"%f\" old_width=\"%f\"/>\n",
            "bottom" . $geneID,
            ($fam_st) * $sz_dif,
            (0) * ($sz_h),
            ($gene_size) * $sz_dif,
            (3) * ($sz_h),
            $back_col[$even],
            ($fam_st_old) * $sz_dif,
            ( $fam_end_old - $fam_st_old ) * $sz_dif
        );

        #Setting the ending value
        $fam_end     = $cnt + $gene_size;
        $fam_end_old = $cnt_old + $gene_size_old;

        #Keeping track of the x-position (cnt),
        $cnt     += $gene_size;
        $cnt_old += $gene_size_old;
        $cur_fam = $f_id;

    }

    #Going through the genes in order to write the gene, fGI names and the lines connecting them
    $cnt = 0;
    for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

        my ( $ID, $num, $type, $f_id, $member, $geneID, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[$j] };
        my $gene_size_old = abs( $end_loc - $st_loc );
        my $gene_size     = $std_gene_x;

        $last_num = $j;

        $fam_end = $cnt + $gene_size;
        $cur_fam = $f_id;

        #removed fgi from top_row, bottom_row
        $top_row .= sprintf(
"<text id=\"%s\" x=\"%f\" y=\"%f\" text-anchor=\"end\" dominant-baseline=\"middle\" transform=\"rotate(270 %f,%f)\" font-size=\"%f\" fam=\"%s\" >%s</text>\n",
            "topText" . $geneID,
            ( $cnt + $gene_size / 2 ) * $sz_dif,
            (1.75) * ($sz_h),
            ( $cnt + $gene_size / 2 ) * $sz_dif,
            (1.75) * ($sz_h),
            $sz_h * 0.7,
            $cur_fam, $f_id
        );
        $bottom_row .= sprintf(
"<text id=\"%s\" x=\"%f\" y=\"%f\" text-anchor=\"start\" dominant-baseline=\"middle\" transform=\"rotate(-90 %f,%f)\" font-size=\"%f\" fam=\"%s\"  onclick=\"writeFasta(\'%s\')\" >%s</text>\n",
            "botText" . $geneID,
            ( $cnt + $gene_size / 2 ) * $sz_dif,
            (3) * ($sz_h),
            ( $cnt + $gene_size / 2 ) * $sz_dif,
            (3) * ($sz_h),
            $sz_h * 0.7,
            $cur_fam, $geneID, $geneID
        );
        $after .= sprintf(
"<line id=\"%s\" x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"darkgray\"  stroke-dasharray=\"0.25,0.75\"/>\n",
            "line" . $geneID,
            ( $cnt + $gene_size / 2 ) * $sz_dif,
            ($i) * ($sz_h),
            ( $cnt + $gene_size / 2 ) * $sz_dif,
            ($i2) * ($sz_h)
        );
        $cnt += $gene_size;

    }

    #String for all the infromation for the left group (bounding 5' gene), right group (bounding 3' gene),
    my $left_grp  = "";
    my $right_grp = "";
    my $r_cnt     = 0;

    #Adding the annotation text for the right side buttons
    $right_grp .= sprintf(
        "<text id=\"%s\" x=\"%f\" y=\"%f\" font-size=\"%f\" transform=\"rotate(-90 %f %f)\">%s</text>",
        "HighlightLeg",
        (0.5) * $sz_h,
        ($dif_i) * $sz_h,
        14,
        (0.5) * $sz_h,
        ($dif_i) * $sz_h, "Highlight"
    );
    $right_grp .= sprintf(
        "<text id=\"%s\" x=\"%f\" y=\"%f\" font-size=\"%f\" transform=\"rotate(-90 %f %f)\">%s</text>",
        "SelectLeg",
        (1.25) * $sz_h,
        ($dif_i) * $sz_h,
        14,
        (1.25) * $sz_h,
        ($dif_i) * $sz_h, "Select"
    );
    $right_grp .= sprintf(
        "<text id=\"%s\" x=\"%f\" y=\"%f\" font-size=\"%f\" transform=\"rotate(-90 %f %f)\">%s</text>",
        "ShowGenomeLeg",
        ( $sz_dif * $max_end + $sz_h * 5.25 ),
        ($dif_i) * $sz_h,
        14,
        ( $sz_dif * $max_end + $sz_h * 5.25 ),
        ($dif_i) * $sz_h,
        "Show Genomes"
    );

    #Going through each of the fGIs to add the different buttons & bounding genes for each of the different gene paths
    for ( my $i2 = 0 ; $i2 < scalar( @{ $data->{order} } ) ; $i2++ ) {

        my $fam_st   = 0;
        my $fam_end  = 0;
        my $cur_fam  = "";
        my $last_num = 0;
        my $i        = $i2 + $dif_i;

        #circle to highlight the row
        $right_grp .= sprintf(
            "<svg id=\"%s\" y=\"%f\" x=\"%f\" height=\"%f\" numRow=\"%f\" ison=\"Off\">\n",
            "All" . $data->{order}->[$i2],
            ($i) * $sz_h,
            0, $sz_h, $i2
        );
        $right_grp .= sprintf(
"<circle id=\"%s\" onclick=\"highlightRow(\'%s\')\" cy=\"%f\" cx=\"%f\" r=\"%f\" fill=\"white\" stroke=\"black\" ison=\"Off\"/>\n",
            "circle" . $data->{order}->[$i2],
            $data->{order}->[$i2],
            (.5) * $sz_h,
            0.25 * $sz_h,
            0.25 * $sz_h
        );

        #Botton to show the genomes or remove genomes below: expands
        $right_grp .= sprintf(
            "<svg id=\"%s\" x_loc=\"%f\" onclick=\"showGenomes(\'%s\')\" y=\"%f\" x=\"%f\">\n",
            "show" . $data->{order}->[$i2],
            $i2, $data->{order}->[$i2],
            0, ( $sz_dif * $max_end + $sz_h * 4.5 )
        );
        $right_grp .= sprintf(
            "<rect y=\"%f\" x=\"%f\" height=\"%f\" width=\"%f\" fill=\"white\" stroke=\"black\"/>\n",
            $sz_h * .25,
            $sz_h * .25,
            $sz_h * .5,
            $sz_h * .5
        );
        $right_grp .= sprintf(
            "<line id =\"%s\" y1=\"%f\" x1=\"%f\" y2=\"%f\" x2=\"%f\" stroke=\"black\"/>\n",
            "showText" . $data->{order}->[$i2],
            $sz_h * .35,
            $sz_h * .5, $sz_h * .65,
            $sz_h * .5
        );
        $right_grp .= sprintf(
            "<line id =\"%s\" y1=\"%f\" x1=\"%f\" y2=\"%f\" x2=\"%f\" stroke=\"black\"/>\n",
            "keepShowText" . $data->{order}->[$i2],
            $sz_h * .5, $sz_h * .35,
            $sz_h * .5, $sz_h * .65
        );

#$right_grp .= sprintf("<text id=\"%s\" x=\"%f\" y=\"%f\" font-size=\"%f\" text-anchor=\"middle\" alignment-baseline=\"central\">+</text>\n", "showText".$data->{order}->[$i2],$sz_h * .5,$sz_h * .5, $sz_h * .65);
        $right_grp .= sprintf("</svg>");

        #Checkbox to turn rows on or off. Used for buttons on bottom right
        $right_grp .= sprintf(
            "<svg id=\"%s\" x_loc=\"%f\" onclick=\"isFastaOn(\'%s\', \'Full\')\" y=\"%f\" x=\"%f\" ison=\"Off\">\n",
            "genomeBox" . $data->{order}->[$i2],
            $i2, $data->{order}->[$i2],
            0, (0.5) * $sz_h
        );
        $right_grp .= sprintf(
            "<rect id=\"%s\" y=\"%f\" x=\"%f\" height=\"%f\" width=\"%f\" fill=\"white\" stroke=\"black\"/>\n",
            "genomeTextBox" . $data->{order}->[$i2],
            $sz_h * .25,
            $sz_h * .25,
            $sz_h * .5, $sz_h * .5
        );

#$right_grp .= sprintf("<text id=\"%s\" x=\"%f\" y=\"%f\" font-size=\"%f\" text-anchor=\"middle\" alignment-baseline=\"central\"> </text>\n", "genomeTextBox".$data->{order}->[$i2],$sz_h * .5,$sz_h * .5, $sz_h * .65);

        $right_grp .= sprintf("</svg>");

        #
        $right_grp .= sprintf(
"<svg id=\"%s\" onclick=\"writeFasta(\'%s\',)\" y=\"%f\" x=\"%f\" height=\"%f\" position=\"absolute\" ison=\"Off\">\n",
            "show" . $i2,
            $data->{st_id}->{ $data->{order}->[$i2] },
            0, 2 * $sz_h, $sz_h
        );
        $right_grp .= sprintf(
            "<text x=\"%f\" y=\"%f\" text-anchor=\"start\" dominant-baseline=\"middle\" font-size=\"%f\" >%s</text>\n",
            $sz_dif * ($max_end) + 1,
            (.5) * ($sz_h),
            $sz_h * 0.8,
            $data->{cnt}->{ $data->{order}->[$i2] }
        );

        #Adding the highlight in the row
        $after .= sprintf(
            "<rect id=\"%s\" num=\"%d\" x=\"%s\" y=\"%s\" height=\"%s\" width=\"%s\" fill=\"%s\" color=\"%s\"/>",
            "highlight" . $data->{order}->[$i2],
            $i2, 0, $i2 * ($sz_h),
            $sz_h, $tot_w, "none", "none"
        );

        $left_grp .= sprintf(
"<svg id=\"%s\" onclick=\"writeFasta(\'%s\')\" visibility=\"visible\"  y=\"%f\" x=\"%f\" height=\"%f\" position=\"absolute\">\n",
            "rowLeft" . $i2,
            $data->{end_id}->{ $data->{order}->[$i2] },
            $i * $sz_h, 0, $sz_h
        );
        my $dif;    #dif is i
        $i = 0;

        #Drawing the 5' bounding the gene
        if ( $data->{st_id}->{ $data->{order}->[$i2] } ) {

            my $st_x = 0;
            $dif = 10;    #The size of the pointing arrow size
            if ( $dif > $sz_dif * $data->{st_sz}->{ $data->{order}->[$i2] } ) {

                $dif = $sz_dif * $data->{st_sz}->{ $data->{order}->[$i2] };

            }

            #is the gene forward or negative facing?
            if ( $data->{st_dir}->{ $data->{order}->[$i2] } == 1 ) {

                $left_grp .= sprintf(
"<path st=\"%f\" end=\"%f\" id=\"%s\" d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" stroke=\"BLACK\" fill=\"%s\" />\n",
                    $data->{st_sz}->{ $data->{order}->[$i2] },
                    $max_st,
                    $data->{order}->[$i2],
                    $sz_dif * ($st_x),
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($max_st) - $dif,
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($max_st),
                    ( $i + .45 ) * ($sz_h),
                    $sz_dif * ($max_st) - $dif,
                    ( $i + .1 ) * ($sz_h),
                    $sz_dif * ($st_x),
                    ( $i + .1 ) * ($sz_h),
                    "#" . $data->{st_col}->{ $data->{order}->[$i2] }
                );

            } else {

                $left_grp .= sprintf(
"<path st=\"%f\" end=\"%f\" id=\"%s\" d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" stroke=\"BLACK\" fill=\"%s\" />\n",
                    $data->{st_sz}->{ $data->{order}->[$i2] },
                    $max_st,
                    $data->{order}->[$i2],
                    $sz_dif * ($st_x),
                    ( $i + .45 ) * ($sz_h),
                    $sz_dif * ($st_x) + $dif,
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($max_st),
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($max_st),
                    ( $i + .1 ) * ($sz_h),
                    $sz_dif * ($st_x) + $dif,
                    ( $i + .1 ) * ($sz_h),
                    "#" . $data->{st_col}->{ $data->{order}->[$i2] }
                );

            }
            $left_grp .= sprintf(
                "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-size=\"%f\">%s</text>\n",
                $sz_dif * ( $st_x + $max_st - $dif ) / 2,
                (.5) * ($sz_h),
                $sz_h * 0.5,
                $data->{st_id}->{ $data->{order}->[$i2] }
            );

        } else    #if no such gene, draw an empty box
        {

            my $st_x = 0;
            $left_grp .= sprintf(
"<rect x=\"%f\" y=\"%f\" height=\"%f\" width=\"%f\" stroke=\"BLACK\" fill=\"white\" stroke-dasharray=\"1,2\"/>\n",
                $sz_dif * ($st_x),
                (.1) * ($sz_h),
                (.8) * ($sz_h),
                $sz_dif * ( $max_st - $st_x ),
                "#" . $data->{st_col}->{ $data->{order}->[$i2] }
            );

            $left_grp .= sprintf(
                "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-size=\"%f\" >%s</text>\n",
                $sz_dif * ( $st_x + $max_st ) / 2,
                (.5) * ($sz_h),
                $sz_h * 0.5, "Break"
            );

        }
        $left_grp .= "</svg>\n";

        #Same for the 3' bounding gene as the 5'
        if ( $data->{end_id}->{ $data->{order}->[$i2] } ) {

            my $st_x = 0;
            $tot = $max_end;
            my $dif = 10;
            if ( $dif > $sz_dif * $data->{end_sz}->{ $data->{order}->[$i2] } ) {

                $dif = $sz_dif * $data->{end_sz}->{ $data->{order}->[$i2] };

            }
            if ( $data->{end_dir}->{ $data->{order}->[$i2] } == 1 ) {

                $right_grp .= sprintf(
                    "<path d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" stroke=\"BLACK\" fill=\"%s\" />\n",
                    $sz_dif * ($st_x),
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($tot) - $dif,
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($tot),
                    ( $i + .45 ) * ($sz_h),
                    $sz_dif * ($tot) - $dif,
                    ( $i + .1 ) * ($sz_h),
                    $sz_dif * ($st_x),
                    ( $i + .1 ) * ($sz_h),
                    "#" . $data->{end_col}->{ $data->{order}->[$i2] }
                );

            } else {

                $right_grp .= sprintf(
                    "<path d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" stroke=\"BLACK\" fill=\"%s\" />\n",
                    $sz_dif * ($st_x),
                    ( $i + .45 ) * ($sz_h),
                    $sz_dif * ($st_x) + $dif,
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($tot),
                    ( $i + .9 ) * ($sz_h),
                    $sz_dif * ($tot),
                    ( $i + .1 ) * ($sz_h),
                    $sz_dif * ($st_x) + $dif,
                    ( $i + .1 ) * ($sz_h),
                    "#" . $data->{end_col}->{ $data->{order}->[$i2] }
                );

            }
            $right_grp .= sprintf(
                "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-size=\"%f\" >%s</text>\n",
                $sz_dif * ( $st_x + $max_end + $dif ) / 2,
                ( $i + .5 ) * ($sz_h),
                $sz_h * 0.5,
                $data->{end_id}->{ $data->{order}->[$i2] }
            );

        } else {

            my $st_x = 0;
            $right_grp .= sprintf(
"<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" stroke=\"BLACK\" fill=\"white\" stroke-dasharray=\"1,2\" num=\"%s\"/>\n",
                $sz_dif * ($st_x),
                ( $i + .1 ) * ($sz_h),
                $sz_dif * ( $max_end - $st_x ),
                (.8) * ($sz_h), $r_cnt++
            );
            $right_grp .= sprintf(
                "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-size=\"%f\">%s</text>\n",
                $sz_dif * ( $st_x + $max_end + $dif ) / 2,
                ( $i + .5 ) * ($sz_h),
                $sz_h * 0.5, "Break"
            );

        }
        $right_grp .= "</svg></svg>\n";

        $cnt = 0;    #keeping track of the gene size for each order (is this kept?)

        for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

            my ( $ID, $num, $type, $f_id, $member, $geneID, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[$j] };

            my $gene_size = abs( $end_loc - $st_loc );
            if ( $member->{ $data->{order}->[$i2] } ) {

                if ( $cur_fam && $cur_fam ne $f_id ) {

                    $cur_fam = $f_id;
                    $fam_st  = $cnt;

                }
                if ( !$cur_fam ) {

                    $cur_fam = $f_id;
                    $fam_st  = $cnt;

                }
                $last_num = $j;
                $fam_end  = $cnt + $gene_size;

            }
            $cnt += $gene_size;

        }

    }

    #Writing all the SVG for the left and right sides of the main window into their respective divisions
    $out .= sprintf(
"<div id=\"leftDiv\" wd=\"%f\" position=\"fixed\"><svg id=\"leftSVG\" version=\"1.2\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"  y=\"0\" x=\"0\" width=\"%f\" height=\"%f\">\n",
        $sz_dif * $max_st,
        $sz_dif * $max_st,
        $sz_h * ($height) + 2
    );
    $out .= "<g id=\"left_border\">" . $left_grp . "</g></svg></div>";
    $out .= sprintf(
"<div id=\"rightDiv\" wd=\"%f\" position=\"fixed\"><svg id=\"rightSVG\" version=\"1.2\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" y=\"0\" x=\"0\" width=\"%f\" height=\"%f\">\n",
        $sz_dif * $max_end + 135,
        $sz_dif * $max_end + 135,
        $sz_h * ($height) + 2
    );
    $out .= "<g id=\"right_border\">" . $right_grp . "</g></svg></div>";

    #writing the main image backgrounds to the SVG
    $out .= $after;

    ###################################################################################
    #Step (6) Adding in the gene images/arrows to the main window
    ###################################################################################

    $cnt = 0;
    my @rows;    #array with the svg for each of the rows (ie a gene : below I initializes
    for ( my $i2 = 0 ; $i2 < scalar( @{ $data->{order} } ) ; $i2++ ) {

        $rows[$i2] = "";

    }
    foreach my $list ( @{ $data->{cords} } ) {

        my ( $ID, $num, $type, $f_id, $member, $geneID, $c, $st_loc, $end_loc, $info ) = @$list;

        my $ht            = ($num) / $GENOME_NUM;
        my $hits          = 0;                            #number of gene paths in which the gene is found... not used
        my $gene_size_old = abs( $end_loc - $st_loc );    #old style gene size
        my $gene_size     = $std_gene_x;
        my $cnt_old       = 0;
        my $dif           = 0;

#drawing an empty rectangle to surround the whole vertical gene area, mainly to allow the title to show on mouse over even when gene is not there/drawn
        $out .= sprintf(
            "<a xlink:title=\"%s\"><rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"none\" stroke=\"none\"/></a>\n",
            $type,
            $cnt * $sz_dif,
            ($dif_i) * ($sz_h),
            $gene_size * $sz_dif,
            ( $dif_i + scalar( @{ $data->{order} } ) ) * $sz_h
        );

        #Going through all the gene paths to draw the arrow-> writes them to rows array
        for ( my $i2 = 0 ; $i2 < scalar( @{ $data->{order} } ) ; $i2++ ) {

            my $i = 0;

            #If the gene is in this gene-path
            if ( $member->{ $data->{order}->[$i2] } ) {

                $hits++;
                $dif = 10;    #again the size of the gene arrow pointer size
                if ( $dif > $gene_size * $sz_dif ) { $dif = $gene_size * $sz_dif; }
                my $cnt2 = $cnt;
                $cnt = 0;     #Not sure why I am resetting $cnt... (x-axis location)
                              #Drawing the gene arrows in svgs with gene specific mouseovers
                if ( $st_loc < $end_loc ) {

                    $rows[$i2] .= sprintf(
"<a xlink:title=\"%s\"><svg id=\"%s\" x=\"%f\"><path d=\"M%f %f L%f %f L%f %f L%f %f L%f %f L%f %f L%f %f Z\" id=\"%s\" fgi=\"%s\" type=\"%s\" fill=\"%s\" stroke=\"BLACK\" onclick=\"writeFasta(\'%s\')\"/></svg></a>\n",
                        $type,
                        ( $i2 . "_" . $geneID ),
                        $cnt2 * $sz_dif,
                        $cnt * $sz_dif,
                        ( $i + 0.5 ) * ($sz_h) - $sz_h * .15,
                        $cnt * $sz_dif + $gene_size * $sz_dif - $dif,
                        ( $i + 0.5 ) * ($sz_h) - $sz_h * .15,
                        $cnt * $sz_dif + $gene_size * $sz_dif - $dif,
                        ( $i + 0.5 ) * ($sz_h) - $sz_h * .45,
                        $cnt * $sz_dif + $gene_size * $sz_dif,
                        ( $i + 0.5 ) * ($sz_h),
                        $cnt * $sz_dif + $gene_size * $sz_dif - $dif,
                        ( $i + 0.5 ) * ($sz_h) + $sz_h * .45,
                        $cnt * $sz_dif + $gene_size * $sz_dif - $dif,
                        ( $i + 0.5 ) * ($sz_h) + $sz_h * .15,
                        $cnt * $sz_dif,
                        ( $i + 0.5 ) * ($sz_h) + $sz_h * .15,
                        $f_id,
                        $geneID . "detail",
                        $data->{order}->[$i2],
                        ( "#" . $c ),
                        $geneID
                      )

                } else {

                    $rows[$i2] .= sprintf(
"<a xlink:title=\"%s\"><svg id=\"%s\" x=\"%f\"><path d=\"M%f %f L%f %f L%f %f L%f %f L%f %f L%f %f L%f %f Z\" id=\"%s\" fgi=\"%s\" type=\"%s\" fill=\"%s\" stroke=\"BLACK\" onclick=\"writeFasta(\'%s\')\"/></svg></a>\n",
                        $type,
                        ( $i2 . "_" . $geneID ),
                        $cnt2 * $sz_dif,
                        $cnt * $sz_dif + $dif,
                        ( $i + 0.5 ) * ($sz_h) - $sz_h * .45,
                        $cnt * $sz_dif + $dif,
                        ( $i + 0.5 ) * ($sz_h) - $sz_h * .15,
                        $cnt * $sz_dif + $gene_size * $sz_dif,
                        ( $i + 0.5 ) * ($sz_h) - $sz_h * .15,
                        $cnt * $sz_dif + $gene_size * $sz_dif,
                        ( $i + 0.5 ) * ($sz_h) + $sz_h * .15,
                        $cnt * $sz_dif + $dif,
                        ( $i + 0.5 ) * ($sz_h) + $sz_h * .15,
                        $cnt * $sz_dif + $dif,
                        ( $i + 0.5 ) * ($sz_h) + $sz_h * .45,
                        $cnt * $sz_dif,
                        ( $i + 0.5 ) * ($sz_h),
                        $f_id,
                        $geneID . "detail",
                        $data->{order}->[$i2],
                        ( "#" . $c ),
                        $geneID
                    );

                }
                $cnt = $cnt2;

            } else {

                $rows[$i2] .= "";

            }

        }
        $cnt     += $gene_size;        #moving to the right down the x-axis the gene_size
        $cnt_old += $gene_size_old;    #moving down the x-axis old style

    }

    #For all the rows draw the svg with the row
    for ( my $i2 = 0 ; $i2 < scalar( @{ $data->{order} } ) ; $i2++ ) {

        my $i = $i2 + $dif_i;
        $out .= sprintf(
            "<svg id=\"%s\" y=\"%f\" old_y=\"%f\" rowID = \"%s\" height=\"%f\" visibility=\"visible\" cnt=\"%d\">%s</svg>",
            "mainRow" . $i2,
            ($i2) * $sz_h,
            ($i2) * $sz_h,
            $data->{order}->[$i2],
            $sz_h, $data->{cnt}->{ $data->{order}->[$i2] },
            $rows[$i2]
        );

    }
    $out .= "</svg>" . $top_row . "</svg>" . $bottom_row . "</svg>";
    $out .= sprintf("</div></body></html>");
    return ($out);

}

#Makes the fGI region for putting it into the preview window pane
sub make_fgi_region_svg {

    my ( $file, $data, $ont, $out ) = @_;
    my $size   = $data->{end} - $data->{st};       #Getting the width of the region in bp
    my $height = scalar( @{ $data->{order} } );    #Getting the height of the region in # of gene paths
    my $reduct = 100;                              #Factor to reduce the size
    my $sz_sp  = 0.00;
    my $gene_width =
      $size / scalar( @{ $data->{cords} } ) / $reduct;    #scaling the gene size to make sure that they are visible
    my $sz   = $gene_height * 0.25;
    my $sz_h = $sz;

    #resize the gene size if it exceeds the entire windows
    if ( $sz * ( 1.5 + scalar( @{ $data->{cords} } ) ) > 2 * ( $border + $max_radius ) ) {

        $sz = 2 * ( $border + $max_radius ) / ( ( 1.5 + scalar( @{ $data->{cords} } ) ) );

    }

    #initializing the preview svg
    $out .= sprintf( "<svg height=\"%f\" width=\"%f\" id=\"$file\" visibility=\"hidden\" overflow=\"visible\" >\n",
        28, $sz_h * ( 2 * $height ) + 2 );

    #drawing the background window
    $out .= sprintf(
"<rect x=\"0\" y=\"0\" width=\"%f\" height=\"%f\" fill=\"white\" stroke=\"black\" grp_id=\"$file\" onclick=\"closeCL(evt)\"/>\n",
        $sz *   ( 1.5 + scalar( @{ $data->{cords} } ) ),
        $sz_h * ($height) + 2
    );

    #Drawing the background for the fGIs
    my $cnt = 0;

    #Going through the rows of gene paths...
    for ( my $i = 0 ; $i < scalar( @{ $data->{order} } ) ; $i++ ) {

        #Variables used to
        my $fam_st   = ( 1 + $sz_sp ) * $sz - $sz * .475;
        my $fam_end  = 0;
        my $cur_fam  = "";
        my $last_num = 0;

        #adding the number of genomes in each path to the edge
        $out .= sprintf(
            "<text x=\"%f\" y=\"%f\" text-anchor=\"end\" font-size=\"%f\">%s</text>\n",
            $sz * ( scalar( @{ $data->{cords} } ) + 1 ),
            ( $i + 1 ) * ($sz_h),
            $sz_h * 0.98,
            $data->{cnt}->{ $data->{order}->[$i] }
        );

        #Going through the columns of the different genes to draw them as needed
        for ( my $j = 1 ; $j <= scalar( @{ $data->{cords} } ) ; $j++ ) {

            my ( $ID, $num, $type, $f_id, $member, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[ $j - 1 ] };

            if ( $member->{ $data->{order}->[$i] } ) {

                #If there has been a change in family: (a) draw a rectangle in the background to show the family
                if ( $cur_fam && $cur_fam ne $f_id ) {

                    $out .= sprintf(
"<rect fam_st=\"$fam_st\" fam_end=\"$fam_end\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fam=\"%s\" fgi=\"%s\" stroke=\"black\" fill=\"lightgray\"/>\n",
                        $fam_st,
                        ( $i + 0.5 ) * ($sz_h) - $sz_h * .45,
                        $fam_end - $fam_st,
                        ( $sz_h * 0.9 ),
                        $cur_fam, $data->{order}->[$i]
                    );
                    $cur_fam = $f_id;
                    $fam_st  = ( $j + $sz_sp ) * $sz - $sz * .475;

                    #And (b) draw a line in the middle of the row connecting the old fGI region to the new region
                    if ( $last_num + 1 < $j ) {

                        $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" fam=\"%s\" fgi=\"%s\" stroke=\"black\" stroke-dasharray=\"1,1\"/>\n",
                            $fam_end, ( $i + 0.5 ) * ($sz_h), $fam_st,
                            ( $i + 0.5 ) * $sz_h, $cur_fam, $data->{order}->[$i]
                        );

                    }

                }

                #This is initializing the fGI region: and drawing the line from the left edge to the box if needed
                if ( !$cur_fam ) {

                    $cur_fam = $f_id;
                    $fam_st  = ( $j + $sz_sp ) * $sz - $sz * .475;
                    if ( $j > 0 ) {

                        $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" fam=\"%s\" fgi=\"%s\" stroke=\"black\" stroke-dasharray=\"1,1\"/>\n",
                            0, ( $i + 0.5 ) * ($sz_h),
                            $fam_st, ( $i + 0.5 ) * $sz_h,
                            $cur_fam, $data->{order}->[$i]
                        );

                    }

                }
                $last_num = $j;
                $fam_end  = ( $j + $sz_sp ) * $sz + $sz * .475;

            }

        }
        if ( !$cur_fam ) {

            my $j = scalar( @{ $data->{cords} } );
            my ( $ID, $num, $type, $f_id, $member, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[ $j - 1 ] };
            $cur_fam  = $f_id;
            $fam_end  = ( $j + $sz_sp ) * $sz + $sz * .475;
            $last_num = $j;

        }

        #Closing the box
        $out .= sprintf(
"<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fam=\"%s\" fgi=\"%s\" stroke=\"black\" fill=\"lightgray\"/>\n",
            $fam_st,
            ( $i + 0.5 ) * ($sz_h) - $sz_h * .45,
            $fam_end - $fam_st,
            ( $sz_h * 0.9 ),
            $cur_fam, $data->{order}->[$i]
        );

        #Drawing a line from the background box to teh right end if there is a gap
        if ( $last_num + 1 < scalar( @{ $data->{cords} } ) ) {

            $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" fam=\"%s\" fgi=\"%s\" stroke=\"black\" stroke-dasharray=\"1,1\"/>\n",
                $fam_end,
                ( $i + 0.5 ) * ($sz_h),
                ( scalar( @{ $data->{cords} } ) - .5 ) * $sz + $sz * .45,
                ( $i + 0.5 ) * $sz_h,
                $cur_fam, $data->{order}->[$i]
            );

        }

    }

    #Going through the genes to draw the colored ovals
    for ( my $j = 1 ; $j <= scalar( @{ $data->{cords} } ) ; $j++ ) {

        my $list = $data->{cords}->[ $j - 1 ];
        my ( $ID, $num, $type, $f_id, $member, $geneID, $c, $st_loc, $end_loc, $info ) = @$list;
        my $ht   = ($num) / $GENOME_NUM;
        my $hits = 0;

        #going through the gene-paths (rows) and if the gene is in the path- draw it
        for ( my $i = 0 ; $i < scalar( @{ $data->{order} } ) ; $i++ ) {

            if ( $member->{ $data->{order}->[$i] } ) {

                $hits++;
                $out .= sprintf(
"<ellipse cx=\"%f\" cy=\"%f\" rx=\"%f\" ry=\"%f\" fam=\"%s\" id=\"%s\" fgi=\"%s\" stroke=\"black\" fill=\"%s\"/>\n",
                    ( $j + $sz_sp ) * $sz,
                    ( $i + 0.5 ) * ($sz_h),
                    $sz * .45, $sz_h * .45,
                    $f_id,
                    $geneID . "detail",
                    $data->{order}->[$i],
                    ( "#" . $c )
                );

            }

        }
        if ($hits) { $cnt++; }

    }

    #closing the svg and returning the string
    $out .= sprintf("</svg>\n");
    return ($out);

}

#Drawing the circular fGR for the main page
sub svg_draw_plasmid_seg {

    my ( $depth, $d1, $d3, $l1, $l2, $c, $href, $name, $id, $type, $out ) = @_;

    #This is which flag to use in the SVG arc string
    my $cir = 0;
    if ( abs( $d3 - $d1 ) > $seq_len / 2 ) { $cir = 1; }
    my $dif = $seq_len / 250;    #This is the max size for the arrow head
    if ( $dif > abs( $d3 - $d1 ) * 2 ) {

        $dif = abs( $d3 - $d1 ) / 2;

    }                            #If it is more than the half size , use the half the size as the arrow head
    my $c2 = $c;
    if ( $type eq "Gene" ) { $c2 = "none"; }
    $id =~ /CL_(\d+)/;
    my $a = $1;

    #Adding in the plasmid javascript for the gene page information
    my $nJscript .= sprintf( "var currentID =\'%s\';\n\nvar allFastaJSON = \'[", $id );

    foreach my $genomeID ( keys( %{ $seqs{byClust}->{$a} } ) ) {

        $nJscript .= sprintf(
            "{\"seq\":\"%s\", \"id\":\"%s\"},",
            $seqs{byClust}->{$a}->{$genomeID},
            $seqs{byClustId}->{$a}->{$genomeID}
        );

    }
    if ( length($nJscript) > 0 ) {

        chop($nJscript);

    }
    $nJscript .= "]\';\n";

    #Adding the json files
    open( OUTJSON, ">", $out_dir . "/$output_id/json/" . $id . ".allFasta.json" );
    print OUTJSON $nJscript;
    close(OUTJSON);

    #5' to 3' gene directions plot uses arced paths
    if ( $d1 < $d3 ) {

        $out->{$depth} .= sprintf(
"<a xlink:title=\"%s\"><path id=\"%s\" loc=\"%f-%f out of %f\" d=\"M%f %f A%f %f 0 $cir,1 %f %f L%f %f L%f %f A%f %f 0 $cir,0 %f %f L %f %f Z\" fill=\"%s\" fill_old=\"%s\" stroke=\"BLACK\"/ type=\"Plasmid\" onclick=\"writeFasta(\'%s\', \'%s\')\"></a>\n",
            $name . "detail", $id, $d1, $d3, $seq_len, screen_x_trans( $d1, $l1 ), screen_y_trans( $d1, $l1 ),
            $l1,              $l1,
            screen_x_trans( $d3 - $dif, $l1 ), screen_y_trans( $d3 - $dif, $l1 ),
            screen_x_trans( $d3, ( $l1 + $l2 ) / 2 ), screen_y_trans( $d3, ( $l1 + $l2 ) / 2 ),
            screen_x_trans( $d3 - $dif, $l2 ), screen_y_trans( $d3 - $dif, $l2 ),
            $l2, $l2,
            screen_x_trans( $d1, $l2 ), screen_y_trans( $d1, $l2 ),
            screen_x_trans( $d1, $l1 ), screen_y_trans( $d1, $l1 ),
            ($c), $c, $id, $href

        );

    }    #3' to 5' gene directions plot uses arced paths
    else {

        $out->{$depth} .= sprintf(
"<a xlink:title=\"%s\"><path id=\"%s\" d=\"M%f %f A%f %f 0 $cir,0 %f %f L%f %f L%f %f A%f %f 0 $cir,1 %f %f L %f %f Z\" fill=\"%s\" fill_old=\"%s\" stroke=\"BLACK\"/ onclick=\"writeFasta(\'%s\',\'%s\')\"></a>\n",
            $name . "detail", $id, screen_x_trans( $d1, $l1 ), screen_y_trans( $d1, $l1 ),
            $l1,              $l1,
            screen_x_trans( $d3 + $dif, $l1 ), screen_y_trans( $d3 + $dif, $l1 ),
            screen_x_trans( $d3, ( $l1 + $l2 ) / 2 ), screen_y_trans( $d3, ( $l1 + $l2 ) / 2 ),
            screen_x_trans( $d3 + $dif, $l2 ), screen_y_trans( $d3 + $dif, $l2 ),
            $l2, $l2,
            screen_x_trans( $d1, $l2 ), screen_y_trans( $d1, $l2 ),
            screen_x_trans( $d1, $l1 ), screen_y_trans( $d1, $l1 ),
            ($c), $c, $id, $href

        );

    }

    #Returning the svg string with the circular fGR drawn
    return ($out);

}

#Returns string containing the whole html page of the core region page
sub make_core_region_page {

    my ( $file, $data, $out ) = @_;

    #Going through each of the genes and getting the faste information for each one
    for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

        my ( $ID, $num, $type, $f_id, $member, $cl_id, $c, $st_loc, $end_loc, $info ) = @{ $data->{cords}->[$j] };
        $num =~ /_(\d+)/;
        my $a = $1;
        my $nJscript = sprintf( "var currentID =\'%s\';\n\nvar allFastaJSON = \'[", $num );

        foreach my $genomeID ( keys( %{ $seqs{byClust}->{$a} } ) ) {

            $nJscript .= sprintf(
                "{\"seq\":\"%s\", \"id\":\"%s\"},",
                $seqs{byClust}->{$a}->{$genomeID},
                $seqs{byClustId}->{$a}->{$genomeID}
            );

        }
        chop($nJscript);
        $nJscript .= "]\';\n";

        #Writing the fasta javascript to the JSON files
        open( OUTJSON, ">", $out_dir . "/$output_id/json/" . $num . ".allFasta.json" );
        print OUTJSON $nJscript;
        close(OUTJSON);

    }
    my $jscript = "";    #The string with the html javascript
    my $lev     = 0;     #Tree level: 0 means no tree loaded

    #Adding the tree as necessary and setting the level to default (4)
    if ( $tree_level > 0 ) {

        $jscript .= "<script type=\"text/javascript\" src=\'../json/tree.json\'></script>\n";
        $lev = 4;

    }

    #Writing the HTML's javacript
    $jscript .= "<script type=\"text/javascript\">
	//<![CDATA[

	";
    my %terms;    #List of the GO terms and functions

    my $multFasta = sprintf("var multiFastaJason=\'{");
    $jscript .= sprintf("\n\nvar fastaJSON = \'{");    #adding the list of to
    my $size   = $data->{end} - $data->{st};           #size of the core region in bp
    my $reduct = 100;                                  #The amount to reduce the region
    my $out_g .= sprintf(
        "<svg height=\"%f\" width=\"%f\" id=\"%s\" visibility=\"hidden\" overflow=\"visible\" >",
        28,
        $gene_height * 3.5 + 2,
        $file . "geneLen",
    );
    my $x1 = 0;                                        #x-level variable for the legend

    #Going through all the genes in the core region
    for ( my $j = 0 ; $j < scalar( @{ $data->{cords} } ) ; $j++ ) {

        my ( $chr, $ID, $st, $end, $def, $type, $len, $c, $info ) = @{ $data->{cords}->[$j] };
        $ID =~ /_(\d+)/;
        my $a = $1;
        $multFasta .= "$a: {";                         #adding gene id (number) to the multifasta

        my $genomeList    = "";                        #string list of the genomes which contains a core gene
        my $genomeLineLen = 0;                         #x-location of the genome list. Used to add breaks to the genome list
        foreach my $genomeID ( keys( %{ $seqs{byClust}->{$a} } ) ) {

            $genomeLineLen += length($genomeID);
            if ( $genomeLineLen > 80 ) {

                $genomeList .= "<br>";
                $genomeLineLen = length($genomeID);

            }                                          #Adding the break if needed in the line is > 80 characters
            $genomeList .= "$genomeID,";

        }

        #Adding the gene information for the gene page to the Javascript string
        $jscript .= sprintf(
"\"%s\":{\"Cluster\":\"%s\", \"Name\":\"%s\",\"# of Genomes\":\"%s\",\"Functional IDs\":\"%s\",\"Associated Terms\":\"%s\",\"Sequence\":\"%s\", \"Maximum Length\":\"%f\", \"Minimum Length\":\"%f\", \"Mean\":\"%f\", \"Standard Deviation\":\"%f\", \"Genomes\":\"%s\"},",
            $ID,                             $ID,                        $clusters{$a}->{protein_name},
            $clusters{$a}->{num_of_members}, $info->{type_ref},          $info->{oth_ref},
            $seqs{cent}->{$a}->{seq},        $geneLenInfo{$a}->{minLen}, $geneLenInfo{$a}->{maxLen},
            $geneLenInfo{$a}->{mean},        $geneLenInfo{$a}->{sd},     $genomeList
        );

        #Adding term and term info to gene ID
        if ( $info->{oth_ref} ) {

            my @term_list = split ";", $info->{oth_ref};
            foreach my $a (@term_list) {

                if ( $ont{$a} ) { $terms{$a} = $ont{$a}; }

            }

        }

    }
    chop($jscript);
    $jscript .= "}\';\n";

    #Making the graph: if needed, will resize to image:

    #Adding the genes and their term ids to the javascript string
    $jscript .= sprintf("\n\nvar termJSON = \'{");
    my $add_to = "";
    foreach my $a ( keys(%terms) ) { $add_to .= sprintf( "\"%s\":\"%s\",", $a, $terms{$a} ); }
    if ($add_to) { chop($add_to); }
    $jscript .= $add_to . "}\';\n";

    $jscript .= "//]]>
</script><script type=\"text/javascript\" src=\"../scripts/fgi.functions.js\"></script>\n";
    $out_g .= "</svg>";

    #This is the starting string for the html
    my $outn = sprintf(
"<html><body xmlns=\"http://www.w3.org/1999/xhtml\"  xmlns:xlink=\"http://www.w3.org/1999/xlink\" onload=\"draw_tree("
          . ( 0 + $lev )
          . ", \'circular\')\">" );

    my $old_gene_height = $gene_height;
    $gene_height *= 2;
    if ( ( $size / $reduct ) < ( $tree_size * 2 ) ) {

        $reduct = $size / ( 2 * $tree_size );

    }

    #adding the javascript string from the above
    $outn .=
        $jscript
      . "<div style=\"top:0; left:0; position: absolute; height: "
      . ( $gene_height * 4.5 )
      . "; width:"
      . ( $size / $reduct + 30 )
      . "; overflow-x:auto; overflow-y:visible;\"  ht=\""
      . ( $gene_height * 6 )
      . "\" wt=\""
      . ( 2 * $tree_size + 45 )
      . "\" baseProfile=\"tiny\"  transRatio=\"1\" id=\"div$file\" >\n";

#adding inputed string which here is the top core region map that is identical to the preview pane (see below for a decription)

    my $id = $file;

    #Drawing the SVG main background
    $outn .= sprintf(
"<svg height=\"%f\" width=\"%f\" style=\"top:0; left:10; position:absolute;\" class=\"$id\" id=\"$file\" visibility=\"hidden\" transRatio=\"1\" overflow=\"visible\" >",
        $gene_height * 4.5,
        $size / $reduct + 30
    );

    #Adding the background
    $outn .= sprintf(
"<rect x=\"1\" y=\"1\" width=\"%f\" height=\"%f\" fill=\"white\" stroke=\"black\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\" onclick=\"closeCL(evt)\"/>",
        $size / $reduct + 28,
        $gene_height * 4 + 2
    );

    #Adding the top background box
    $outn .= sprintf(
"<rect x=\"%f\" y=\"%f\" width = \"%f\" height=\"%f\" fill=\"lightgray\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
        15,
        $gene_height * 0.125,
        $size / $reduct,
        $gene_height * 1.25
    );

    #Adding the bottom background box
    $outn .= sprintf(
"<rect x=\"%f\" y=\"%f\" width = \"%f\" height=\"%f\" fill=\"lightgray\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
        15,
        $gene_height * 2.375,
        $size / $reduct,
        $gene_height * 1.25
    );

    #Adding the top background box
    $outn .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"black\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
        15,
        $gene_height * 1.5,
        $size / $reduct + 15,
        $gene_height * 1.5
    );

    #Print the ticks
    my $tick_size = 10000;
    for ( my $i = ( int( $data->{st} / $tick_size ) + 1 ) * $tick_size ; $i < $data->{end} ; $i += $tick_size ) {

        #adding the in the top and bottom ticks as well as the BP counts with commas
        $outn .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"black\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.75,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.5
        );
        $outn .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"white\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 0.125,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.375
        );

        $outn .= sprintf(
            "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" font-size=\"14\" transRatio=\"1\" class=\"$id\" >%s</text>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.75 + 8,
            commas($i)
        );

        $outn .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"black\" class=\"$id\" transRatio=\"1\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 2.5,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 2.25
        );
        $outn .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"white\" class=\"$id\" transRatio=\"1\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 3.5,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 2.375
        );

    }

    #Print the genes;
    my $i = 1;
    my %legend;    #hash that matches function name with the color in the legend
    my $sz_dif = 0.05;    #

    foreach my $list ( @{ $data->{cords} } ) {

        my ( $chr, $ID, $st, $end, $def, $type, $len, $c, $info ) = @$list;

        $ID =~ /CL_(\d+)/;
        my $ID2 = $1;     #ID2 is the gene number

        #dif is the x-distance of the arrow point
        my $dif = 10;

        #resets the x-distance of the arrow if it is currently longer than the length
        if ( abs( $end - $st ) / $reduct < 10 ) {

            $dif = abs( $end - $st ) / ( 2 * $reduct );

        }

        #Sets the color is there isn't a function
        if ( !$c ) { $c = "#bbbbbb"; }
        my $new_c = $c;       # color of the gene
        $new_c =~ s/\s//g;    #removing any spaces from the color
        if ( !$new_c ) { $new_c = "#000000"; }    #setting default value to black
        $legend{ $info->{type_ref} } = $new_c;    #setting the color of on the legend
        if ( $st < $end ) {

            #This is the core region gene arrow
            $outn .= sprintf(
"<a xlink:title=\"%s\"><path id=\"%s\" transRatio=\"1\" class=\"$id\" d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"writeFasta(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010; Mean Length: "
                  . $geneLenInfo{$ID2}->{mean}
                  . " &#177; "
                  . $geneLenInfo{$ID2}->{sd},
                $ID . "detail",
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 0.25,
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct - $dif + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct + 15,
                $gene_height * 0.75,
                ( $end - $data->{st} ) / $reduct - $dif + 15,
                $gene_height * 0.25,
                ($c),
                $ID
            );

            #This is the hieght of the core region
            my $ht = ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM;

            #Printing out the bar plot showing the relative number of genomes with the core gene as shown
            $outn .= sprintf(
"<a xlink:title=\"%s\"><rect id=\"%s\" transRatio=\"1\" class=\"$id\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"score_tree(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010;# of Genomes: "
                  . int( $clusters{$ID2}->{num_of_members} ),
                $ID . "detailHT",
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * ( 2.50 + ( ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM ) * 4 ),
                ( $end - $st ) / $reduct,
                $gene_height * ( ( $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM - $CORE_PERCT ) * 4,
                ($c),
                $ID
            );
            my $n_lines = int( $clusters{$ID2}->{num_of_members} / $reduct + 0.5 );

        } else {

            $outn .= sprintf(
"<a xlink:title=\"%s\"><path id=\"%s\" transRatio=\"1\" class=\"$id\" d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"writeFasta(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010;Mean Length: "
                  . $geneLenInfo{$ID2}->{mean}
                  . " &#177; "
                  . $geneLenInfo{$ID2}->{sd},
                $ID . "detail",
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 0.25,
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct + $dif + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct + 15,
                $gene_height * 0.75,
                ( $end - $data->{st} ) / $reduct + $dif + 15,
                $gene_height * 0.25,
                ($c),
                $ID
            );

            my $ht = ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM;

            $outn .= sprintf(
"<a xlink:title=\"%s\"><rect id=\"%s\" x=\"%f\" transRatio=\"1\" class=\"$id\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"score_tree(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010;# of Genomes: "
                  . int( $clusters{$ID2}->{num_of_members} ),
                $ID . "detailHT",
                ( $end - $data->{st} ) / $reduct + 15,
                $gene_height * ( 2.50 + ( ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM ) * 4 ),
                ( $st - $end ) / $reduct,
                $gene_height * ( ( $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM - $CORE_PERCT ) * 4,
                ($c),
                $ID
            );

        }

    }
    $outn .= "</svg></div>\n";

    my $tot_width = $tree_size * 2;
    my $sz_h      = $gene_height * 0.5;
    my @leg_id    = keys(%legend);        #All the functions in the legend
    my $x1n       = 0;                    #location on the x-axis
    my $y_len     = 1;                    #hieght of this particular rows
    my @max_n;                            #array that stores the Maximum number of text lines in each of the horizontal rows
    my @x_len;                            #array that stores the x location of each legend key
    my @y_loc;                            #array that stores the y location of each legend key
    my $x_l = 0;                          #legend text horizontal size. Reset at each number
    my @y_l;                              #array that keeps track of the y location os
    my $cur_y         = 0;                #The current vertical level of the legend
    my $legend_width  = 0;                # the width of the legend on the bottom
    my $legend_height = 0;                #The height of the legend
    my $sz_h2         = $sz_h * 0.5;      #Another gene height that is set at half the height

    #Getting the x and y location on the legend for each of the functions in the fGI
    #Also helping to set up the svg
    my $r1 = 0.25;    #Horizontal size of a letter
    $max_n[ $y_len - 1 ] = 0;    #
    for ( my $i = 0 ; $i < scalar(@leg_id) ; $i++ ) {

        my $y1      = $sz_h * ( $gene_height * 5 ) + 2;                   #The vertical location
        my $max_l   = 0;                                                  #Maximum vertival size in this row
        my $id_name = split_id( $leg_id[$i], $x1n + $sz_h * 0.6, 15 );    #divides an id into different lines as needed
        my $n       = -1;
        while ( $id_name =~ /([^\>\<]+)\<\/tspan\>/g )                    #Getting the number of lines in a in id
        {

            $n++;
            if ( length($1) > $max_l ) { $max_l = length($1); }           # re

        }

        $x_l = $max_l * $sz_h * $r1 + $sz_h * 0.75;
        if ( $x1n + $x_l >
            $tot_width )  #if the legend exceeds the horizontal size, this moves down to lower y-value and resets the x-value
        {

            $x_len[ $y_len - 1 ] = $x1n;    #resetting x-value to the start
            $legend_height += ( 1 + $max_n[ $y_len - 1 ] );            #adding to the legend height
            $cur_y += $sz_h * ( ( $max_n[ $y_len - 1 ] + 1 ) * 1 );    #
            if ( $x1n + $x_l > $legend_width ) { $legend_width = $x1n + $x_l + 4; }    #resetting the legend width as needed
            $y_len++;
            $x1n = 0;

        }
        if (  !$max_n[ $y_len - 1 ]
            || $n >
            $max_n[ $y_len - 1 ] )    #if the number in a line exceeds the current maximium, makes thje max equal this number
        {

            $max_n[ $y_len - 1 ] = $n;

        }
        $x1n += $x_l;                 #Adding to x value
        $y_loc[$i] = $y_len - 1;      #Keeping y as needed
        $y_l[$i]   = $cur_y;

    }

    $x_len[ $y_len - 1 ] = $x1n;
    $legend_height += $cur_y + ( $max_n[ $y_len - 1 ] + 1 );
    if ( $x1n > $legend_width ) { $legend_width = $x1n + 4; }

    #Setting the x location of the legend box
    $x1n = ( $tot_width - $legend_width ) / 2;
    if ( $x1n < 0 ) {

        $x1n          = 0;
        $legend_width = $tot_width * $sz_dif;

    }    #If x1 is greater than the total width, then just set it to zero
         #$legend_width = max(250, $legend_width);
    $outn .= sprintf(
"<div id = \"legendDiv\" ht=\"%f\" wd=\"%f\"  style=\"height:%f; width:%f; position:absolute; top:%f;\" num=\"%f\"><svg id=\"legendSVG\" version=\"1.2\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id=\"$file\" height=\"\%f\" width=\"\%f\" left=\"\%f\">\n",
        ( $legend_height + 1.5 ) * $sz_h, $legend_width, ( $legend_height + 1.5 ) * $sz_h,
        $tot_width,                       $gene_height * 5, scalar(@leg_id),
        ( $legend_height + 1.5 ) * $sz_h, $legend_width,    $x1n
    );
    $x1n = 0;
    $outn .= sprintf("<g id=\"legend\">\n");

    #Adding background rectangle
    $outn .=
      sprintf( "<rect id=\"legRect\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"#eeeeee\" stroke =\"black\"/>",
        $x1n, 0, $legend_width, ( $legend_height + 1.5 ) * $sz_h );

    #Adding the text title
    $outn .= sprintf(
        "<text id=\"LegText\" x=\"%f\" y=\"%f\" text-anchor=\"middle\" font-size=\"%s\">Legend</text>\n",
        $x1n + $legend_width / 2,
        $sz_h * ( $legend_height + 1.25 ) + 2,
        $sz_h * 1.25
    );

    $y_len = 0;

    #Writing all the functional annotation
    for ( my $i = 0 ; $i < scalar(@leg_id) ; $i++ ) {

        my $x_l     = 0;
        my $id_name = split_id( $leg_id[$i], $x1n + $sz_h * 0.6, 15 );    #the functional name as a multiple line
        my $max_l   = 0;
        my $n       = -1;
        while ( $id_name =~ /([^\>\<]+)\<\/tspan\>/g ) {

            $n++;
            if ( length($1) > $max_l ) { $max_l = length($1); }

        }

        $x_l = $max_l * $sz_h * $r1 + $sz_h * 0.75;                       #getting the width of the legend text

        if ( $y_loc[$i] != $y_len ) { $y_len++; $x1n = ( $tot_width - $x_len[$y_len] ) / 2; }
        my $y1 = 0;
        if ( $y_l[$i] ) { $y1 = $y_l[$i]; }
        $id_name = split_id( $leg_id[$i], $sz_h * 0.6, 15 );              #resetting the id name by removing the
               #Writing out the legend functional text and colored square in an SVG
        $outn .= sprintf( "<svg id=\"%s\"  x=\"%f\" y=\"%f\">", "TextID" . $i, $x1n, $y1 );
        $outn .= sprintf(
            "<rect  href=\"%s\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" />\n",
            $i, 0,
            $sz_h * ( $max_n[$y_len] ) / 5,
            $sz_h * 0.5,
            $sz_h / 3.0,
            $legend{ $leg_id[$i] }
        );
        $outn .= sprintf(
"<text href=\"%s\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"auto\" font-size=\"%f\" dominant-baseline=\"hanging\" text-anchor=\"start\" fill=\"black\">%s</text>\n",
            $i,
            $sz_h * 0.5,
            $sz_h * ( $max_n[$y_len] - $n ) / 5,
            $border * 1.15,
            $sz_h * .5, $id_name
        );
        $outn .= "</svg>";
        $x1n += $x_l;

    }
    $out .= "</g>\n";

    #Adding in the resize buttons

    #Adding in the svg_g
    $outn .= $out_g;
    $outn .= sprintf( "</div><div id=\"sortDiv\" style=\"top:%fpx; left:1px; position: relative;\">",
        ( $legend_height + 1.5 ) * $sz_h + $gene_height * 5 );
    $outn .= sprintf(
"<button style=\"left:1; top:5; position: absolute; width: %f; height: %f; background-color:#4444ff; \" onclick=\"transformGenes(\'$file\', \'95\')\"  id=\"decreaseButton\"> - </button>",
        20, $gene_height / 3 );
    $outn .=
      sprintf( "<strong style=\"left:26; top:5; position: absolute; width: %f; height: %f;\" > Change View Zoom </strong>",
        140, $gene_height / 3 );
    $outn .= sprintf(
"<button style=\"left:171; top:5; position: absolute; width: %f; height: %f; background-color:#4444ff; \" onclick=\"transformGenes(\'$file\', \'105\')\"  id=\"increaseButton\"> + </button>",
        20, $gene_height / 3 );
    $outn .= "</div>";

    #Adding in the tree and the tree functions when possible
    if ( $tree_level > 0 ) {

        #This is the disk images for the saving the tree as SVG/PNG
        my $trans          = sprintf( "translate(%f, %f) scale(%f, %f)", 0, 30, 0.045, 0.045 );
        my $action         = "onclick=\"saveSVGtree(\'svg\', \'$file\')\"";
        my $disk_image_svg = $disk_image_orig;
        $disk_image_svg =~ s/TRANS/$trans/g;
        $disk_image_svg =~ s/ACTION/$action/g;
        $disk_image_svg =~ s/TEXT/SVG/g;
        $disk_image_svg =~ s/FS/250/g;

        $trans = sprintf( "translate(%f, %f) scale(%f, %f)", 60, 30, 0.045, 0.045 );
        $action = "onclick=\"saveSVGtree(\'png\', \'$file\')\"";
        my $disk_image_png = $disk_image_orig;
        $disk_image_png =~ s/TRANS/$trans/g;
        $disk_image_png =~ s/ACTION/$action/g;
        $disk_image_png =~ s/TEXT/PNG/g;
        $disk_image_png =~ s/FS/250/g;

        #This is making the table that serves to divide the tree screen using $tree_size
        $outn .= sprintf(
"\n<div style=\"top:%fpx; left:0; position: relative;\"><table ><tr height=\"%fpx\"><td style=\"vertical-align:top;\" rowspan=\"3\"><div><svg version=\"1.2\" overflow=\"visible\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id=\"tree_svg\" height=\"%fpx\" width=\"%fpx\"></svg></div></td><td style=\"vertical-align:top;\"><div><h3>Tree Functions</h3>Outer Ring Tree Level:<button onclick=\"changeTree(\'-\')\">-</button><b id = \"LevelCount\">$lev</b><button onclick=\"changeTree(\'+\')\">+</button><form>Phylogeny Style:<input type=\"radio\" name=\"type\" value=\"Circular\" checked=\"true\" onclick=\"changeType(event)\"/>Circular<input type=\"radio\" name=\"type\" value=\"Linear\" onclick=\"changeType(event)\"/>Linear</form><svg \"id=\"ConstantLegendSVG\" version=\"1.2\" overflow=\"visible\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" id=\"tree_save_svg\" x=\"0\" height=\"%fpx\" width=\"%fpx\"><text x=\"2\" y=\"14\">Node Colors:</text><circle cx=\"95\" r=\"4\" cy=\"8\" fill=\"black\"/><text x=\"102\" y=\"14\">Present</text><circle cx=\"165\" r=\"4\" cy=\"8\" fill=\"LightGray\"/><text x=\"172\" y=\"14\">Missing</text>"
              . $disk_image_svg
              . $disk_image_png
              . "</svg></div></td></tr><tr height=\"%fpx\" ><td style=\"vertical-align:top;\"><div id=\"gene_list\" style=\"height:%fpx; overflow-y: auto;\"></div></td></tr><tr height=\"%fpx\"><td style=\"vertical-align:bottom;\"><div id=\"percentage_image\"></div></td></tr><tr height=\"%fpx\" ><td><table><tr><td height = \"%fpx\" width=\"%fpx\" id=\"Select Metatype\" style=\"vertical-align:top;\"><h3>Metadata Type</h3>",
            ( $legend_height + 1.5 ) * $sz_h + $gene_height * 6,
            $tree_size * 0.2,
            $tree_size,
            $tree_size,
            $tree_size * 0.1,
            $tree_size * 0.4,
            $tree_size * 0.4,
            $tree_size * 0.4,
            $tree_size * 0.15,
            $tree_size * 0.25,
            $tree_size * 0.2,
            ( $tree_size * 0.4 )
        );

        #Getting a list of the metadata associated with the genomes (added by user)
        my @k = keys(%meta_types);

        #Adding a checkbox to for each metadata for each
        for ( my $i = 0 ; $i < scalar(@k) ; $i++ ) {

            my $a = $k[$i];
            $outn .= sprintf(
                "<input type=\"checkbox\" id=\"select$a\" value=\"$a\" onchange=\"change_metadata_status(event)\" />$a<br />"
            );

        }

        #Adding the legend
        $outn .= sprintf(
"</td><td id = \"Legend\" style=\"vertical-align:top;\"><h3>Legend</h3><svg id=\"LegendSVG\" version=\"1.2\" overflow=\"visible\" baseProfile=\"tiny\" xmlns=\"http://www.w3.org/2000/svg\" height=\"%f\" width=\"%f\"></svg></td></tr></table>",
            $tree_size * 0.15,
            $tree_size * 0.6
        );
        $outn .= "</td><td style=\"vertical-align:top;\"><h3>Percentages</h3><div id =\"percDiv\"></div></td></tr></table>";
        $outn .= "</div>";

    }

    #closing the html
    $outn .= "</body></html>\n";

    #returning the string with the html
    $gene_height = $old_gene_height;
    return $outn;

}

#Making the window pane preview of the core region. This is also used in the core region full page
sub make_core_region_svg {

    my ( $file, $id, $data, $ont, $out ) = @_;
    my $size   = $data->{end} - $data->{st};    #Getting the size in bp
    my $reduct = 100;                           #The proportion to reduce: a constant

    #Drawing the SVG main background
    $out .= sprintf(
"<svg height=\"%f\" width=\"%f\" class=\"$id\" id=\"$file\" visibility=\"hidden\" transRatio=\"1\" overflow=\"visible\" >",
        $gene_height * 4.5,
        $size / $reduct + 30
    );

    #Adding the background
    $out .= sprintf(
"<rect x=\"1\" y=\"1\" width=\"%f\" height=\"%f\" fill=\"white\" stroke=\"black\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\" onclick=\"closeCL(evt)\"/>",
        $size / $reduct + 28,
        $gene_height * 4 + 2
    );

    #Adding the top background box
    $out .= sprintf(
"<rect x=\"%f\" y=\"%f\" width = \"%f\" height=\"%f\" fill=\"lightgray\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
        15,
        $gene_height * 0.125,
        $size / $reduct,
        $gene_height * 1.25
    );

    #Adding the bottom background box
    $out .= sprintf(
"<rect x=\"%f\" y=\"%f\" width = \"%f\" height=\"%f\" fill=\"lightgray\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
        15,
        $gene_height * 2.375,
        $size / $reduct,
        $gene_height * 1.25
    );

    #Adding the top background box
    $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"black\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
        15,
        $gene_height * 1.5,
        $size / $reduct + 15,
        $gene_height * 1.5
    );

    #Print the ticks
    my $tick_size = 10000;
    for ( my $i = ( int( $data->{st} / $tick_size ) + 1 ) * $tick_size ; $i < $data->{end} ; $i += $tick_size ) {

        #adding the in the top and bottom ticks as well as the BP counts with commas
        $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"black\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.75,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.5
        );
        $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"white\" transRatio=\"1\" class=\"$id\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 0.125,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.375
        );

        $out .= sprintf(
            "<text x=\"%f\" y=\"%f\" text-anchor=\"middle\" font-size=\"8\" transRatio=\"1\" class=\"$id\" >%s</text>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 1.75 + 8,
            commas($i)
        );

        $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"black\" class=\"$id\" transRatio=\"1\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 2.5,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 2.25
        );
        $out .= sprintf(
"<line x1=\"%f\" y1=\"%f\" x2 = \"%f\" y2=\"%f\" stroke=\"white\" class=\"$id\" transRatio=\"1\" grp_id=\"$file\"/>\n",
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 3.5,
            ( $i - $data->{st} ) / $reduct + 15,
            $gene_height * 2.375
        );

    }

    #Print the genes;
    my $i = 1;
    foreach my $list ( @{ $data->{cords} } ) {

        my ( $chr, $ID, $st, $end, $def, $type, $len, $c, $info ) = @$list;

        $ID =~ /CL_(\d+)/;
        my $ID2 = $1;    #ID2 is the gene number

        #dif is the x-distance of the arrow point
        my $dif = 10;

        #resets the x-distance of the arrow if it is currently longer than the length
        if ( abs( $end - $st ) / $reduct < 10 ) {

            $dif = abs( $end - $st ) / ( 2 * $reduct );

        }

        #Sets the color is there isn't a function
        if ( !$c ) { $c = "#bbbbbb"; }
        if ( $st < $end ) {

            #This is the core region gene arrow
            $out .= sprintf(
"<a xlink:title=\"%s\"><path id=\"%s\" transRatio=\"1\" class=\"$id\" d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"writeFasta(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010; Mean Length: "
                  . $geneLenInfo{$ID2}->{mean}
                  . " &#177; "
                  . $geneLenInfo{$ID2}->{sd},
                $ID . "detail",
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 0.25,
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct - $dif + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct + 15,
                $gene_height * 0.75,
                ( $end - $data->{st} ) / $reduct - $dif + 15,
                $gene_height * 0.25,
                ($c),
                $ID
            );

            #This is the hieght of the core region
            my $ht = ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM;

            #Printing out the bar plot showing the relative number of genomes with the core gene as shown
            $out .= sprintf(
"<a xlink:title=\"%s\"><rect id=\"%s\" transRatio=\"1\" class=\"$id\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"score_tree(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010;# of Genomes: "
                  . int( $clusters{$ID2}->{num_of_members} ),
                $ID . "detailHT",
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * ( 2.50 + ( ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM ) * 4 ),
                ( $end - $st ) / $reduct,
                $gene_height * ( ( $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM - $CORE_PERCT ) * 4,
                ($c),
                $ID
            );
            my $n_lines = int( $clusters{$ID2}->{num_of_members} / $reduct + 0.5 );

        } else {

            $out .= sprintf(
"<a xlink:title=\"%s\"><path id=\"%s\" transRatio=\"1\" class=\"$id\" d=\"M%f %f L%f %f L%f %f L%f %f L%f %f Z\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"writeFasta(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010;Mean Length: "
                  . $geneLenInfo{$ID2}->{mean}
                  . " &#177; "
                  . $geneLenInfo{$ID2}->{sd},
                $ID . "detail",
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 0.25,
                ( $st - $data->{st} ) / $reduct + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct + $dif + 15,
                $gene_height * 1.25,
                ( $end - $data->{st} ) / $reduct + 15,
                $gene_height * 0.75,
                ( $end - $data->{st} ) / $reduct + $dif + 15,
                $gene_height * 0.25,
                ($c),
                $ID
            );

            my $ht = ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM;

            $out .= sprintf(
"<a xlink:title=\"%s\"><rect id=\"%s\" x=\"%f\" transRatio=\"1\" class=\"$id\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"%s\" stroke=\"BLACK\" stroke-width=\"1\" grp_id=\"$file\" onclick=\"score_tree(\'%s\')\"/></a>\n",
                "Name: "
                  . $clusters{$ID2}->{protein_name}
                  . "&#010;# of Genomes: "
                  . int( $clusters{$ID2}->{num_of_members} ),
                $ID . "detailHT",
                ( $end - $data->{st} ) / $reduct + 15,
                $gene_height * ( 2.50 + ( ( $GENOME_NUM - $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM ) * 4 ),
                ( $st - $end ) / $reduct,
                $gene_height * ( ( $clusters{$ID2}->{num_of_members} ) / $GENOME_NUM - $CORE_PERCT ) * 4,
                ($c),
                $ID
            );

        }

    }
    $out .= "</svg>\n";
    return ($out);

}

#Drawing an SVG tick using svg_spoke onto the $out string
sub svg_tick {

    my ( $x, $l, $out ) = @_;
    my ( $h1, $h2 );

    $h1 = convert_layer($l);

    $h2 = $h1 + $tick_height;

    $out = svg_spoke( 100, $x, $h1, $h2, "#000000", 0, 0, 0, $out );
    return ($out);

}

#Drawing an arc dsegment
sub svg_draw_arc_seg {

    my ( $depth, $d1, $d3, $l1, $l2, $c, $href, $name, $id, $type, $out ) = @_;
    $out->{$depth} .= "<g>\n";
    $d1 += $cur_st;
    $d3 += $cur_st;
    if ($href) { $out->{$depth} .= "<a xlink:title=\"$name\">\n"; }
    my $d2 = ( ( $d3 - $d1 ) / 2 ) + $d1;
    if ($name) { $out->{$depth} .= "<title>$name<\/title>\n"; }

    #onmouseout=\"turnPlotOff(evt)\"
    my $c2 = $c;
    if ( $type && $type eq "Gene" ) { $c2 = "none"; }
    my $func  = "";
    my $id_in = $id;
    if ( $depth == 1 ) {

        $func = "onclick=\"loadSVG(evt)\" onmouseout=\"turnPlotOff(evt)\" onmouseover=\"showPlot(evt)\"";

    }
    if ( $depth == 3 ) {

        $func  = "onclick=\"writeFasta(\'$name\', \'Plasmid$id\')\"";
        $id_in = $href;

    }
    if ( ( $d3 - $d1 ) / $seq_len < 0.5 ) {

        $out->{$depth} .= sprintf(
            "<path d=\"M %f %f\nA %f,%f 0 0,1 %f,%f\n",
            screen_x_trans( $d1, $l1 ),
            screen_y_trans( $d1, $l1 ),
            $l1, $l1,
            screen_x_trans( $d3, $l1 ),
            screen_y_trans( $d3, $l1 )
        );
        $out->{$depth} .= sprintf(
"L %f %f\nA %f,%f 0 0,0 %f,%f\nL %f %f\" id=\"$id_in\" stroke = \"%s\" fill = \"%s\" type = \"$type\" fill_old = \"%s\" href=\"$href\" $func angle=\"%s\"\/>\n",
            screen_x_trans( $d3, $l2 ),
            screen_y_trans( $d3, $l2 ),
            $l2,
            $l2,
            screen_x_trans( $d1, $l2 ),
            screen_y_trans( $d1, $l2 ),
            screen_x_trans( $d1, $l1 ),
            screen_y_trans( $d1, $l1 ),
            $c2,
            $c,
            $c,
            360 * ( ( $d1 + $d2 ) / ( 2 * $seq_len ) )
        );

    } else {

        $out->{$depth} .= sprintf(
            "<path d=\"M %f %f\nA %f,%f 0 1,1 %f,%f\n",
            screen_x_trans( $d1, $l1 ),
            screen_y_trans( $d1, $l1 ),
            $l1, $l1,
            screen_x_trans( $d3, $l1 ),
            screen_y_trans( $d3, $l1 )
        );
        $out->{$depth} .= sprintf(
"L %f %f\nA %f,%f 0 1,0 %f,%f\nL %f %f\" id=\"$id_in\" stroke = \"%s\" fill = \"%s\" type = \"$type\" fill_old = \"%s\" href=\"$href\" $func angle=\"%s\"\/>\n",
            screen_x_trans( $d3, $l2 ),
            screen_y_trans( $d3, $l2 ),
            $l2,
            $l2,
            screen_x_trans( $d1, $l2 ),
            screen_y_trans( $d1, $l2 ),
            screen_x_trans( $d1, $l1 ),
            screen_y_trans( $d1, $l1 ),
            $c2,
            $c,
            $c,
            360 * ( ( $d1 + $d2 ) / ( 2 * $seq_len ) )
        );

    }
    if ($href) { $out->{$depth} .= "</a>\n"; }
    $out->{$depth} .= "<\/g>\n";
    return ($out);

}

sub svg_draw_arc_full {

    my ( $depth, $l1, $l2, $c, $href, $name, $out ) = @_;
    if ($href) { $out->{$depth} .= "<a xlink:href=\"$href\" xlink:title=\"$name\">\n"; }
    $out->{$depth} .= sprintf(
        "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"%s\" fill=\"none\" stroke-width=\"%f\"/>",
        get_center_x(), get_center_y(), ( ( $l2 + $l1 ) / 2 ),
        $c, abs( $l2 - $l1 )
    );
    if ($href) { $out->{$depth} .= "</a>\n"; }
    return $out;

}

sub place_svg_ticks {

    my ( $increment, $size, $out ) = @_;
    my ( $h, $i, $j, $d, $justification );

    # this puts a tick mark at intervals around the circle

    #Writing ticks on level 100
    for ( $i = 0 ; $i < $size ; $i += $increment ) {

        $out->{100} = svg_tick( $i, 0, $out->{100} );

        $d = 360 * ( $i / $seq_len );
        $justification = 0;
        if ( $d > 180 ) {

            $justification = 2;

        }

        if ( $d > 120 && $d < 240 ) {

            $h = convert_layer(0) + $tick_height + $text_height;

        } else {

            $h = convert_layer(0) + $tick_height + ( $text_height / 2 );

        }

        $j = $i;
        $j = 1 if ( $j == 0 );

        # THIS IS WHERE LABELS ARE PLACED:
        ### due to the different locations around the pie, modification is
        ### required
        my $degrees = ( 360 * ( $j / $seq_len ) );
        my $lngth = ( length( commas($j) ) * 5 ) + 5;
        my ( $mod_x, $mod_y );
        if ( ( $degrees > 0 ) && ( $degrees <= 90 ) ) {

            $mod_x = 0;
            $mod_y = 0;    #-5;  ## height of one of the character

        }
        if ( ( $degrees > 90 ) && ( $degrees <= 180 ) ) {

            $mod_x = 0;
            $mod_y = 0;

        }
        if ( ( $degrees > 180 ) && ( $degrees <= 270 ) ) {

            $mod_x = 0;    #-$lngth;
            $mod_y = 0;

        }
        if ( ( $degrees > 270 ) && ( $degrees <= 360 ) ) {

            $mod_x = 0;    #-$lngth;
            $mod_y = 0;    #-5; ### again the height

        }
        my $anchor = "start";

        my $dn = $d + 270;
        if ( $d > 180 ) { $anchor = "end"; $dn -= 180; }
        if ( $dn > 360 ) { $dn -= 360; }
        $out->{100} .= sprintf(
"<text x=\"%f\" y=\"%f\" fill=\"black\" transform=\"rotate(%f,%f,%f)\" text-anchor=\"$anchor\" style=\"dominant-baseline: central;\">%s<\/text>\n",
            screen_x_trans( $j, $h ) + $mod_x,
            screen_y_trans( $j, $h ) + $mod_y,
            $dn,
            screen_x_trans( $j, $h ) + $mod_x,
            screen_y_trans( $j, $h ) + $mod_y,

            commas($j)
        );

    }
    return ($out);

}

sub place_svg_ticks_arc {

    my ( $increment, $size, $out ) = @_;
    my ( $h, $i, $j, $d, $justification, $mod_x, $mod_y );

    # this puts a tick mark at intervals around the circle

    #Writing ticks on level 100
    for ( $i = $cur_st ; $i < $size + $cur_st ; $i += $increment ) {

        $out->{100} = svg_tick( $i, 0, $out->{100} );

        $d = 360 * ( $i / $seq_len );
        $justification = 0;
        if ( $d > 180 ) {

            $justification = 2;

        }

        if ( $d > 120 && $d < 240 ) {

            $h = convert_layer(0) + $tick_height + $text_height;

        } else {

            $h = convert_layer(0) + $tick_height + ( $text_height / 2 );

        }

        $j = $i - $cur_st;
        $j = 1 if ( $j == 0 );

        # THIS IS WHERE LABELS ARE PLACED:
        ### due to the different locations around the pie, modification is
        ### required
        my $degrees = ( 360 * ( $j / $seq_len ) );
        my $lngth = ( length( commas($j) ) * 5 ) + 5;
        if ( ( $degrees > 0 ) && ( $degrees <= 90 ) ) {

            $mod_x = 0;
            $mod_y = 0;    #-5;  ## height of one of the character

        }
        if ( ( $degrees > 90 ) && ( $degrees <= 180 ) ) {

            $mod_x = 0;
            $mod_y = 0;

        }
        if ( ( $degrees > 180 ) && ( $degrees <= 270 ) ) {

            $mod_x = 0;    #-$lngth;
            $mod_y = 0;

        }
        if ( ( $degrees > 270 ) && ( $degrees <= 360 ) ) {

            $mod_x = 0;    #-$lngth;
            $mod_y = 0;    #-5; ### again the height

        }
        my $anchor = "start";

        my $dn = $d + 270;
        if ( $d > 180 ) { $anchor = "end"; $dn -= 180; }
        if ( $dn > 360 ) { $dn -= 360; }

        #Writing out the tick labels, including with commas
        $out->{100} .= sprintf(
"<text x=\"%f\" y=\"%f\" fill=\"black\" transform=\"rotate(%f,%f,%f)\" text-anchor=\"$anchor\" style=\"dominant-baseline: central;\">%s<\/text>\n",
            screen_x_trans( $i, $h ) + $mod_x,
            screen_y_trans( $i, $h ) + $mod_y,
            $dn,
            screen_x_trans( $i, $h ) + $mod_x,
            screen_y_trans( $i, $h ) + $mod_y,

            commas($j)
        );

    }
    return ($out);

}

sub place_svg_ticks_multi {

    my ( $increment, $tick, $out ) = @_;
    my ( $h, $i, $j, $d, $justification, $mod_x, $mod_y );

    my @chr = keys(%$tick);
    $i = 0;

    #Writing ticks on level 100
    for ( my $i1 = 0 ; $i1 < scalar(@chr) ; $i1++ ) {

        my $h1 = convert_layer(0) + $tick_height;

        my $h2 = convert_layer(1) - $tick_height;

        $out->{100} = svg_spoke( 0, $tick->{ $chr[$i1] }->{st}, $h1, $h2, "#000000", 0, 0, 0, $out->{100} );
        if ( $increment > 0 ) {

            for ( my $q = $increment ; $q < $tick->{ $chr[$i1] }->{sz} ; $q += $increment ) {

                $i = $tick->{ $chr[$i1] }->{st} + $q;
                $out = svg_tick( $i, 0, $out );

                $d = 360 * ( $i / $seq_len );

                $justification = 0;
                if ( $d > 180 ) {

                    $justification = 2;

                }

                if ( $d > 120 && $d < 240 ) {

                    $h = convert_layer(0) + $tick_height + $text_height;

                } else {

                    $h = convert_layer(0) + $tick_height + ( $text_height / 2 );

                }

                $j = $i;
                $j = 1 if ( $j == 0 );

                # THIS IS WHERE LABELS ARE PLACED:
                ### due to the different locations around the pie, modification is
                ### required
                my $degrees = ( 360 * ( $j / $seq_len ) );
                my $lngth = ( length( commas($q) ) * 5 ) + 5;
                if ( ( $degrees > 0 ) && ( $degrees <= 90 ) ) {

                    $mod_x = 0;
                    $mod_y = 0;    #-5;  ## height of one of the character

                }
                if ( ( $degrees > 90 ) && ( $degrees <= 180 ) ) {

                    $mod_x = 0;
                    $mod_y = 0;

                }
                if ( ( $degrees > 180 ) && ( $degrees <= 270 ) ) {

                    $mod_x = 0;    #-$lngth;
                    $mod_y = 0;

                }
                if ( ( $degrees > 270 ) && ( $degrees <= 360 ) ) {

                    $mod_x = 0;    #-$lngth;
                    $mod_y = 0;    #-5; ### again the height

                }
                my $anchor = "start";

                my $dn = $d + 270;
                if ( $d > 180 ) { $anchor = "end"; $dn -= 180; }
                if ( $dn > 360 ) { $dn -= 360; }
                $out->{100} .= sprintf(
"<text x=\"%f\" y=\"%f\" fill=\"black\" transform=\"rotate(%f,%f,%f)\" text-anchor=\"$anchor\" style=\"dominant-baseline: central;\">%s<\/text>\n",
                    screen_x_trans( $j, $h ) + $mod_x,
                    screen_y_trans( $j, $h ) + $mod_y,
                    $dn,
                    screen_x_trans( $j, $h ) + $mod_x,
                    screen_y_trans( $j, $h ) + $mod_y,
                    commas($q)
                );

            }

        }

    }
    return ($out);

}

sub svg_spoke {

    my ( $depth, $d, $l1, $l2, $c, $width, $href, $name, $out ) = @_;

    # width is optional

    if ($href) { $out .= "<a xlink:href=\"h:\\$href\" xlink:title=\"$name\">\n"; }

    $out .= sprintf(
        "<path d=\"M %d %d\nL %d %d\" stroke=\"black\" angle=\"%f\"\/>\n",
        screen_x_trans( $d, $l1 ),
        screen_y_trans( $d, $l1 ),
        screen_x_trans( $d, $l2 ),
        screen_y_trans( $d, $l2 ),
        360 * ( $d / $seq_len )
    );
    if ($href) { $out .= "</a>\n"; }
    return ($out);

}

sub convert_layer {

    my ($l) = @_;
    return ( $max_radius - ( $layer_height * $l ) );

}

sub screen_x_trans {

    my ( $p, $r ) = @_;
    my ( $d, $x );

    $d = 360 * ( $p / $seq_len );
    return ( ( $r * sin( deg2rad($d) ) ) + $max_radius + $border );

}

sub screen_y_trans {

    my ( $p, $l ) = @_;
    my ($d);

    $d = 360 * ( $p / $seq_len );
    return ( ( $l * cos( deg2rad($d) ) * -1 ) + $max_radius + $border );

}

sub screen_x_trans_ratio {

    my ( $p, $r, $sz ) = @_;
    my ( $d, $x );

    #my $ratio = $sz / (2.0 * $max_radius + 2.0 * $border);
    $d = 360 * ( $p / $seq_len );
    return ( ( $r * sin( deg2rad($d) ) + $sz / 2 ) );

}

sub screen_y_trans_ratio {

    my ( $p, $l, $sz ) = @_;
    my ($d);

    $d = 360 * ( $p / $seq_len );

    #my $ratio = $sz / (2.0 * $max_radius + 2.0 * $border);

    return ( ( $l * cos( deg2rad($d) ) * -1 ) + $sz / 2 );

}

sub rad2deg {

    my ($r) = @_;

    return ( ( ( 180 / $PI ) * ($r) ) % 360 );

}

sub deg2rad {

    my ($d) = @_;

    return ( ( $d / 180 ) * $PI );

}

sub get_center_x {

    return ( $max_radius + $border );

}

sub get_center_y {

    return ( $max_radius + $border );

}

sub commas {    # Adds commas to large numbers to clean up printing the coordingates around the circle and for any other use as needed 

    my ($_) = @_;
    1 while s/(.*\d)(\d\d\d)/$1,$2/;

    $_;

}

#Returns a string to initaite
sub svg_init_image {

    my ( $x, $y, $ret ) = @_;
    $ret->{beg} = sprintf(
"<svg version=\"1.2\" baseProfile=\"tiny\" width=\"$x\" height=\"$y\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id = \"allsvg\"><defs id=\"defs\"></defs>\n"
    );
    $ret->{end} = "</svg>\n";
    return ( %{$ret} );

}

#Returns the bigger of two values
sub max {

    my ( $x, $y ) = @_;
    return ( $x >= $y ) ? $x : $y;

}

#Returns the smaller of two values
sub min {

    my ( $x, $y ) = @_;
    return ( $x < $y ) ? $x : $y;

}

#Returning an html tspan where the provided words can be in multiple rows
sub split_id {

    my $word = $_[0];
    my $x    = $_[1];
    my $len  = $_[2];
    my @n    = split " ", $word;

    #$out is the string
    my $out = "<tspan x=\"$x\" dy=\"0\" href=\"$word\">";
    my $cur = 0;
    foreach my $a (@n) {

        #Trims the word by changing "and" to a plus sign
        if ( $a eq "and" ) { $a = "+"; }
        if ( $a eq "\&" )  { $a = "+"; }

        #makes a new line if it gets beyond the length
        if ( length($a) + $cur > $len ) {

            $out .= "</tspan><tspan href=\"$word\" x=\"$x\" dy=\"1em\">$a";
            $cur = length($a);

        } else {

            #Otherwise just add it to the next line
            $out .= " $a";
            $cur += length($a);

        }

    }

    #Returning the tspan string with a
    return ( $out . "</tspan>" );

}

sub make_main_javascript {

    my $SVGHEIGHT_BIG = $SVGHEIGHT * $DPI;
    my $SVGWIDTH_BIG  = $SVGWIDTH * $DPI;

    #This is the disk image for the saving both for the fGI page and the core tree page
    my $trans = sprintf( "translate(%f, %f) scale(%f, %f)", 0, 30, 0.045, 0.045 );
    my $action = "onclick=\\\"saveSVGtree(\\\'svg\\\')\\\"";

    #disk image with SVG written on it
    my $disk_image_svg = $disk_image_js;
    $disk_image_svg =~ s/TRANS/$trans/g;
    $disk_image_svg =~ s/ACTION/$action/g;
    $disk_image_svg =~ s/TEXT/SVG/g;
    $disk_image_svg =~ s/FS/250/g;

    $trans = sprintf( "translate(%f, %f) scale(%f, %f)", 60, 30, 0.045, 0.045 );
    $action = "onclick=\\\"saveSVGtree(\\\'png\\\')\\\"";

    #disk image with PNG written on it
    my $disk_image_png = $disk_image_js;
    $disk_image_png =~ s/TRANS/$trans/g;
    $disk_image_png =~ s/ACTION/$action/g;
    $disk_image_png =~ s/TEXT/PNG/g;
    $disk_image_png =~ s/FS/250/g;


	
    open( OUT, ">", $out_dir . "/$output_id/scripts/main.javascript.js" );

    print OUT "

//###################################################################
//Javascript functions used by the PanACEA Main Page
//Written by Thomas Clarke
//###################################################################


// Variable telling whether the program is activly computing
var isWorking = \"no\";

// Whether the legend is clickable- 
var LegendOn = \"yes\";

// Array of table ids that show whether the gene is on
var good;

//Array of 
var reg2gene;

var treeJSON = null; //Tree JSON is either null if there is no phylogeny or an array if there is

var searchIsOn = null;

/*Function that initiates the main page by:
	(a) loading JSON information passed in the html file into variables
	(b) starting the loading screen and waiting screen, but hiding the later
	(c) setting the main page to the 1st chromosome
*/

function init()
{
	//Setting the main page as visible
	var divMain = document.getElementById(\"mainDiv\");
	divMain.style.visibility = \"visible\";
	console.log(divMain.style.visibility);
	good = JSON.parse(goodJSON);
	reg2gene = JSON.parse(region_2_gene);
	document.body.style.cursor = \"default\";

	//After all is loaded: hide the loading and waiting pages
	var divWait = document.getElementById(\"Loading\");
	divWait.style.visibility=\"hidden\";
	//var wait = document.getElementById(\"WaitingDiv\");
	//wait.style.visibility = \"hidden\";

	//Setting all the animation duration on loading/waiting page to zero
	var load = document.getElementsByClassName(\'twirl\');
	for (var i = 0; i < load.length; i++)
	{
		load[i].style.animationDuration = \"0s\";
	}
	var load = document.getElementsByClassName(\'loader\');
	for (var i = 0; i < load.length; i++)
	{
		load[i].style.animationDuration = \"0s\";
	}
	var load = document.getElementsByClassName(\'shifting\');
	for (var i = 0; i < load.length; i++)
	{
		load[i].style.animationDuration = \"0s\";
	}

	//Setting the main page to the first chromosome
	setSVG(1);
}

function showSearch(type)
{
	var svg_all = document.getElementById(\"allsvg\");
	if (!svg_all.hasAttribute(\"searchVal\"))
	{
		if (searchIsOn == null)
		{
			var glass = document.getElementById(\"magnifyGlass\");
			glass.setAttribute(\"fill\", \"yellow\");
			var searchBox = document.getElementById(\"searchBox\");
			searchBox.style.visibility = \"visible\";
			searchIsOn = \"On\";
			var searchRect = document.getElementById(\"searchRect\");
			searchRect.style.fill=\"DarkGray\";
		}
		else
		{
			if (type === \"icon\")
			{
				var glass = document.getElementById(\"magnifyGlass\");
				glass.setAttribute(\"fill\", \"yellow\");
				var searchBox = document.getElementById(\"searchBox\");
				searchBox.style.visibility = \"hidden\";
				searchIsOn = null;
				var searchRect = document.getElementById(\"searchRect\");

				searchRect.style.fill=\"LightGray\";
			}
		}
	}
}

function searchInit(evt)
{
	evt.preventDefault();
	runType(\"runSearch\");
}

function searchOn()
{
	var glass = document.getElementById(\"magnifyGlass\");
	glass.setAttribute(\"fill\", \"yellow\");
	var searchBox = document.getElementById(\"searchBox\");
	searchBox.style.visibility = \"visible\";
	
}


function searchOff()
{
	var glass = document.getElementById(\"magnifyGlass\");
	var searchBox = document.getElementById(\"searchBox\");
	if (searchIsOn == null)
	{
		glass.setAttribute(\"fill\", \"black\");
		searchBox.style.visibility =  \"hidden\";
	}
}

/*Shows the appropriate background during a mouse move
Used to make sure that the user knows the program is working
*/
function mouseMove()
{
	var divMain = document.getElementById(\"mainDiv\");

	if (isWorking == \"no\")
	{
		divMain.style.opacity = \"1.0\";

	}
	if (isWorking == \"yes\")
	{
		divMain.style.opacity = \"0.5\";
	}
}

/*When the program is \"working\" the following is turned on:
	(a) puts the curser to wait
	(b) places an arc moving around circle in the middle\
	(c) dims the rest of the page
	(d) sets the global variable isWorking to yes
*/
function Working()
{
	var divMain = document.getElementById(\"mainDiv\");
	var divWait = document.getElementById(\"Waiting\");
	divWait.style.visibility=\"visible\";
	console.log(\"Is Working..\");
	//var wait = document.getElementById(\"WaitingDiv\");
	//wait.style.visibility = \"visible\";

	//Setting cursor values
	divMain.style.pointerEvents = \"none\";
	document.body.style.cursor = \"wait\";

	//setting the animation time for the waiting screen to 2 seconds
	var load = document.getElementsByClassName(\'loader\');
	for (var i = 0; i < load.length; i++)
	{
		load[i].style.animationDuration = \"2s\";
	}

	//setting the animation time for the waiting screen to 3 seconds
	var load = document.getElementsByClassName(\'shifting\');
	for (var i = 0; i < load.length; i++)
	{
		load[i].style.animationDuration = \"3s\";
	}
	isWorking=\"yes\";
}

/*When the program has stopped working, the following is turned off:
	(a) the curser goes back to normal
	(b) the moving arc is hidden
	(c) the page is no longer dimmed
	(d) the isWorking is set back to no
*/
function Done()
{
	var divWait = document.getElementById(\"Waiting\");
	divWait.style.visibility=\"hidden\";
	//var wait = document.getElementById(\"WaitingDiv\");
	//wait.style.visibility = \"hidden\";
	var divMain = document.getElementById(\"mainDiv\");
	divMain.style.opacity = \"1.0\";
	divMain.style.pointerEvents = \"auto\";
console.log(\"Is Done..\");
	document.body.style.cursor = \"default\";
	var load = document.getElementsByClassName(\'loader\');
	for (var i = 0; i < load.length; i++)
	{
		load[i].style.animationDuration = \"0s\";
	}
	var load = document.getElementsByClassName(\'shifting\');
	for (var i = 0; i < load.length; i++)
	{
		load[i].style.animationDuration = \"0s\";
	}
	isWorking=\"no\";

}

/*Wrapper function to handle any number of clicks, mouse overs, etc
	Used to prevent multiple actions being taken at once
	Passed variables: type = name of function to use,
						evt = event variable,
						id= string with name of id
*/
function runType(type, evt, id)
{
	//only will work if something else is going on
	if (isWorking == \"no\")
	{
		if (type != \"GeneType\" || LegendOn == \"yes\")
		{
		//function to start the working
		Working();
		//setTimeout is used to make sure that working is turned on
		setTimeout(function(){

		//
		
		if (type ==\"resetSearch\")
		{
			resetSearch();
		}
		if (type ==\"runSearch\")
		{
			console.log(\"Running search...\")
			searchTable();
		}
		if (type == \"showRowPlot\")
		{
			showRowPlot(evt);
		}

		//
		if (type == \"GeneType\")
		{
			if (LegendOn == \"yes\")
			{
				selectGeneType(evt);
			}
		}

		//
		if (type == \"TableStatus\")
		{
			changeTableStatus(evt);
		}

		//
		if (type == \"TableText\")
		{
			saveTableTxt(id);
		}

		//
		if (type == \"TableButton\")
		{
			selectTableButton(evt, id);
		}

		//
		if (type == \"SelectMenu\")
		{
			selectMenu(evt);
		}

		//
		if (type == \"ChangeMenu\")
		{
			changeMenu(evt);
		}


		Done();},0);
		}
	}
}

//Wrapper to run turnSampleOff
function turnPlotOff(evt) {
	var target = evt.target;
	turnSampleOff(target);
}

//function to open file in the ID href from the html
function loadSVG(evt)
{
	var target =evt.target;
	if (target.hasAttribute(\"href\"))
	{
		var addr = target.getAttribute(\"href\");
		var location = window.location.href;
		var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);

		window.open(curPath + addr);
	}
}

/*Function that sets the main/central image to the chromosome with the ID
passed to it (num). Used by init() and by running the side menu
*/
function setSVG(num)
{
	var svg_all = document.getElementById(\"allsvg\");

	//Checks to see if the main page is already turned on
	if (svg_all.hasAttribute(\"coreNum\"))
	{
		var oldNum = svg_all.getAttribute(\"coreNum\");

		//Only do something if the selected number isn't already on the main page
		if (oldNum != num)
		{

			var newSVG = document.getElementById(\"svgOut\"+num);
			var newPreviewSVG = document.getElementById(\"svgPreview\"+num);
			var newArcSVG = document.getElementById(\"svgArc\"+num);

			if (newSVG == null )
			{
				//put up an error popup if there SVG couldn't be found
				alert(\"Cannot find \" +num);
			}
			else
			{
				//Hides the previously shown chromosome
				var oldSVG = document.getElementById(\"svgOut\"+oldNum);
				oldSVG.setAttribute(\"visibility\", \"hidden\");
				var oldPreviewSVG = document.getElementById(\"svgPreview\"+oldNum);
				oldPreviewSVG.setAttribute(\"visibility\", \"hidden\");
				var oldArcSVG = document.getElementById(\"svgArc\"+oldNum);
				oldArcSVG.setAttribute(\"visibility\", \"hidden\");

				//Turns the new ones on
				newSVG.setAttribute(\"visibility\", \"visible\");
				newPreviewSVG.setAttribute(\"visibility\", \"visible\");
				newArcSVG.setAttribute(\"visibility\", \"visible\");

				//Sets the Core Number
				svg_all.setAttribute(\"coreNum\", num);
				if (svg_all.hasAttribute(\"tableType\"))
				{
					makeTable(svg_all.getAttribute(\"tableType\"));
				}
			}
		}
	}
	else
	{

		//Initializes the main page with a chromosome
		var newSVG = document.getElementById(\"svgOut\"+num);
		var newArcSVG = document.getElementById(\"svgArc\"+num);
		var newPreviewSVG = document.getElementById(\"svgPreview\"+num);
		if (newSVG == null)
		{
			alert(\"Cannot find \" +num);
		}
		else
		{
			newSVG.setAttribute(\"visibility\", \"visible\");
			newArcSVG.setAttribute(\"visibility\", \"visible\");

			newPreviewSVG.setAttribute(\"visibility\", \"visible\");

			//Turn on the coreNum and Table variables
			svg_all.setAttribute(\"coreNum\", num);
			if (svg_all.hasAttribute(\"tableType\"))
			{
				makeTable(svg_all.getAttribute(\"tableType\"));
			}
		}
	}
}

function selectGeneChr()
{
	var svg_all = document.getElementById(\"allsvg\");
	var table = null;
	var tableHead = JSON.parse(tableHeadStr);
	var tableWidth = JSON.parse(tableWidthStr);
	var tableInfo = JSON.parse(tableInfoStr);
	var kys = Object.keys(good);
	
}


/*Function to allow the user to open or close sections of the left hand side menu
*/
function selectMenu(evt)
{
	//target1 => Core Menu with thumbnails is the moved out version
	//target1 => Core Menu with thumbnails is the moved in version

	var target1 = document.getElementById(\"coreMenuOut\");
	var target2 = document.getElementById(\"coreMenuIn\");

	//Getting the main SVG
	var svg_all = document.getElementById(\"allsvg\");

	//Seeing if the menu is \"on\"
	if (svg_all.hasAttribute(\"menuOn\"))
	{
		//Turning off the menu (ie moving it to a closed state)
		target1.style.left = \"0px\";
		target2.style.width = \"0px\";

		//Menu on flag is turned off
		svg_all.removeAttribute(\"menuOn\");
		var thumb = target2.childNodes;

		//Hiding all the thumbnails
		for (var i = 0; i < thumb.length; i++)
		{
			if (thumb[i].tagName==\"svg\")
			{
				thumb[i].setAttribute(\"visibility\", \"hidden\");
			}
		}

		//Get the blank SVG (empty space used to buffer the menu)
		var thumb = document.getElementById(\"blankSVG\");
		if (thumb != null)
		{
			thumb.setAttribute(\"visibility\", \"hidden\");
		}
	}
	else
	{
		//set the thumb size to the user set values
		target2.style.width = $thumb_size_w + \"px\";
		target1.style.left = $thumb_size_w + \"px\";

		//Turn menu on variable on
		svg_all.setAttribute(\"menuOn\", \"1\");
		var thumb = target2.childNodes;

		//Turning all the thumb
		for (var i = 0; i < thumb.length; i++)
		{
			if (thumb[i].tagName==\"svg\")
			{

				thumb[i].setAttribute(\"visibility\", \"visible\");
			}
		}

		//Same for blank buffer site
		var thumb = document.getElementById(\"blankSVG\");
		if (thumb != null)
		{
			thumb.setAttribute(\"visibility\", \"visible\");
		}
	}
}

/*Changes the thumbnail menu by showing or hiding the groups of thumbnail sketches of the pan-chromosomes
*/
function changeMenu(evt)
{
	var target = evt.target;
	var typeID = target.getAttribute(\"type\");
	var target2 = document.getElementById(\"coreMenuIn\");
	var thumb = target2.childNodes;

	//Going through the the sketches and either hiding them (height = 0) and showing them
	for (var i = 0; i < thumb.length; i++)
	{
		if (thumb[i].tagName==\"svg\" && thumb[i].getAttribute(\"type\") == typeID)
		{
			if (target.getAttribute(\"Showing\") == \"on\")
			{
				thumb[i].style.height=\"0\";
			}
			else
			{
				thumb[i].style.height=\"$thumb_size_h\";
			}
		}
	}
	//Changing the menu header to show that it is off (having a plus sign that shows that it can be expanded)
	if (target.getAttribute(\"Showing\") == \"on\")
	{
		target.setAttribute(\"Showing\", \"off\");
		target.innerHTML = target.id + \" (+)\";
	}
	else
	{
		//Changing the menu header to show that it is on (having a minus sign that shows that in can be contracted)
		target.setAttribute(\"Showing\", \"on\");
		target.innerHTML = target.id + \" (-)\";
	}
}

/*
Turns the table row off and all the core regions
*/
function turnSampleOff(target)
{
	var id_name = target.getAttribute(\"href\");
	var center = document.getElementById(id_name);
	var svg_all = document.getElementById(\"allsvg\");

	if (center != null)
	{
		center.setAttribute(\"visibility\", \"hidden\");
	}
	target.setAttribute(\"fill\", target.getAttribute(\"fill_old\"));
	svg_all.removeAttribute(\"core_button\");
	svg_all.removeAttribute(\"core_set\");
}

function restoreLegend()
{
	var Leg1 = document.getElementById(\"Legend1\");
	var rects = Leg1.querySelectorAll(\".LegRect\");
	var texts = Leg1.querySelectorAll(\".LegText\");
	var cols = Leg1.querySelectorAll(\".LegCol\");
	for (var i =0; i < rects.length; i++)
	{
		rects[i].setAttribute(\"opacity\", \"1.0\");
	}
	for (var i =0; i < texts.length; i++)
	{
		texts[i].setAttribute(\"opacity\", \"1.0\");
	}
	for (var i =0; i < cols.length; i++)
	{
		cols[i].setAttribute(\"opacity\", \"1.0\");
	}
}

function turnRed(evt, newCol)
{
	var target = evt.target;
	if (newCol != \"old\")
	{
		target.setAttribute(\"oldcol\", target.style.stroke)
		target.style.stroke = newCol;
	}
	else
	{
		target.style.stroke = target.getAttribute(\"oldcol\");
		target.removeAttribute(\"oldcol\");
	}
}

function resetSearch()
{
	var svg_all = document.getElementById(\"allsvg\");
	var table = null;
	var tableInfo = JSON.parse(tableInfoStr);
	var showText = document.getElementById(\"oldSearch\");
	document.getElementById(\"tableDiv\").removeChild(showText);
	document.getElementById(\"searchForm\").style.visibility = \"\";
	document.getElementById(\"searchText\").value = \"\";
	var kys = Object.keys(tableInfo);
	searchBox.style.opacity = \"1\";

	for (var keyName in kys)
	{
		var idType = kys[keyName];
		for (i = 0; i < tableInfo[idType].length; i++)
		{
			good[idType][tableInfo[idType][i][\"href\"]] = 1;
		}
	}
	if (svg_all.hasAttribute(\"tableOn\"))
	{
		makeTable(svg_all.getAttribute(\"tableType\"));
		if (svg_all.hasAttribute(\"expand\"))
		{
			expandTable();
		}
	}
	for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
	{
		if (\"detail\" in tableInfo[\"Gene\"][i])
		{
			var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
			if (gene != null)
			{
				gene.setAttribute(\"fill\", gene.getAttribute(\"fill_old\"));
				gene.setAttribute(\"stroke\", \"none\");
				gene.setAttribute(\"fill-opacity\", \"1.0\");
				gene.setAttribute(\"stroke-opacity\", \"1.0\");
			}
			var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
			if (gene != null)
			{
				gene.setAttribute(\"fill-opacity\", \"1.0\");
				gene.setAttribute(\"stroke-opacity\", \"1.0\");
			}
			var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
			if (gene != null)
			{
				gene.setAttribute(\"fill-opacity\", \"1.0\");
				gene.setAttribute(\"stroke-opacity\", \"1.0\");
			}
		}
	}
	for (var i = 0; i < tableInfo[\"Region\"].length; i++)
	{		
		var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
		good[\"Region\"][tableInfo[\"Region\"][i][\"href\"]] = 1;
		if (gene != null)
		{

			gene.setAttribute(\"fill-opacity\", \"1.0\");
			gene.removeAttribute(\"stroke\");
		}
	}
	svg_all.removeAttribute(\"searchVal\");
	
	showSearch();
}	

function searchTable()
{
	var svg_all = document.getElementById(\"allsvg\");
	var table = null;

	var showText = document.getElementById(\"searchText\").value;
	var searchText = showText.toLowerCase();
	
	
	//Getting the table information from the JSON strings
	var tableHead = JSON.parse(tableHeadStr);
	var tableInfo = JSON.parse(tableInfoStr);
	var term2type = JSON.parse(term_2_type);

	if (svg_all.hasAttribute(\"highlight.row.id\"))
	{
		var turnOff = document.getElementById(svg_all.getAttribute(\"highlight.row.id\"));
		turnOff.style.backgroundColor = target.getAttribute(\"old_bg\");
		var turnOffPart = document.getElementById(turnOff.getAttribute(\"href\"));
		if (turnOffPart != null)
		{
			turnSampleOff(turnOffPart);
		}
		svg_all.removeAttribute(\"highlight.row.id\");
	}
	
	var selectType = \"\";
	var othType = \"\";
	console.log(searchText);
	var totalHits = 0;
	var searchBar = document.getElementById(\"searchForm\");
	var searchIcon = document.getElementById(\"searchIcon\");
	var searchBox = document.getElementById(\"searchBox\");
	var tableDiv = document.getElementById(\"tableDiv\");
	
	var widthVal = parseFloat(searchIcon.style.width)+10;
	
	svg_all.setAttribute(\"searchVal\", searchText);
	
	var showText = document.createElementNS(\"http://www.w3.org/2000/svg\", \"svg\");
	showText.setAttribute(\"id\", \"oldSearch\");
	showText.style.position = \"absolute\";
	showText.style.height = searchIcon.style.height;
	showText.style.width = \"100%\";
	showText.style.top = \"0%\";
	showText.style.left = widthVal +\"px\";
	
	
	showText.setAttribute(\"onclick\", \"runType(\\\"resetSearch\\\")\");
	var rectText = document.createElementNS(\"http://www.w3.org/2000/svg\", \"rect\");
	rectText.setAttribute(\"height\", \"100%\");
	rectText.setAttribute(\"width\", \"100%\");
	rectText.setAttribute(\"left\", \"0%\");
	
	rectText.setAttribute(\"fill\", \"#F5DEB3\");
	
	var textRect = document.createElementNS(\"http://www.w3.org/2000/svg\", \"text\");
	textRect.setAttribute(\"x\", \"10\");
	textRect.setAttribute(\"y\", \"50%\");
	textRect.setAttribute(\"alignment-baseline\", \"middle\");
	textRect.setAttribute(\"font-color\", \"black\"); 
	textRect.innerHTML = searchText;
	
	
	
	
	searchBox.style.opacity = \"0\";
	showText.appendChild(rectText);
	showText.appendChild(textRect);
	tableDiv.appendChild(showText);
	
	var textWidth = textRect.getBoundingClientRect().width;
	var textHeight = textRect.getBoundingClientRect().height;

	console.log(textWidth);
	rectText.setAttribute(\"width\", textWidth+60);
	rectText.setAttribute(\"rx\", \"15\");
	rectText.setAttribute(\"ry\", \"15\");
	
	var textHeight = showText.getBoundingClientRect().height;
	
	var textClose = document.createElementNS(\"http://www.w3.org/2000/svg\", \"path\");
	textClose.setAttribute(\"d\", \"M \" + (textWidth +15) + \",\" + 5 + \" L \" + (textWidth +35) + \",\" + (textHeight-5) + \" M \" + (textWidth +35) + \",\" + 5 + \" L \" + (textWidth +15) + \",\" + (textHeight-5));
	textClose.setAttribute(\"style\", \"stroke: black; stroke-width\: 4;\");
	textClose.setAttribute(\"onmouseover\", \"turnRed(event, \\\"red\\\")\");
	textClose.setAttribute(\"onmouseout\", \"turnRed(event, \\\"old\\\")\");
	showText.style.width = textWidth+70;
	showText.appendChild(textClose);
	
	searchBar.style.visibility = \"hidden\";
	
	var OthTerms = [searchText];
	var kys = Object.keys(tableInfo);
	var searchRect = document.getElementById(\"searchRect\");

	searchRect.style.fill=\"LightGray\";
	for (var keyName in kys)
	{
		var idType = kys[keyName];
		if (idType != \"Gene\" && idType != \"Region\")
		{
			for (var i = 0; i < tableInfo[idType].length; i++)
			{
				var hit_count = 0;
				var ObjKeys = Object.keys(tableInfo[idType][i]);
				for (var j in ObjKeys)
				{
					if (typeof(tableInfo[idType][i][ObjKeys[j]]) === 'string')
					{
						if (tableInfo[idType][i][ObjKeys[j]].toLowerCase().indexOf(searchText) > -1)
						{
							hit_count = hit_count + 1;
						}
					}
				}
				if (hit_count > 0)
				{
					OthTerms.push(tableInfo[idType][i][\"href\"].toLowerCase());
					good[idType][tableInfo[idType][i][\"href\"]] = 1;
					console.log(idType + \" \" + tableInfo[idType][i][\"href\"]);
				}
				else
				{
					good[idType][tableInfo[idType][i][\"href\"]] = 0;
				}
			}
	
		}
	}
	
	for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
	{
		var hit_count =0;
		var ObjKeys = Object.keys(tableInfo[\"Gene\"][i]);

		for (var j in ObjKeys)
		{
			if (typeof(tableInfo[\"Gene\"][i][ObjKeys[j]]) === 'string')
			{
				for (var k in OthTerms)
				{
					if (tableInfo[\"Gene\"][i][ObjKeys[j]].toLowerCase().indexOf(OthTerms[k]) > -1)
					{
						hit_count = hit_count +1;
					}
				}
			}
		}
		if (hit_count > 0)
		{
			totalHits = totalHits + 1;
			good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] = 1;
			
			if (\"detail\" in tableInfo[\"Gene\"][i])
			{
					//Coloring turned on genes
				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
				if (gene != null)
				{
					gene.setAttribute(\"fill\", \"blue\");
					gene.setAttribute(\"stroke\", gene.getAttribute(\"fill_old\"));
				}
			}
			if (\"Region\" in tableInfo[\"Gene\"][i])
			{
				
				var regions = tableInfo[\"Gene\"][i][\"Region\"].split(\",\");
				for (var j in regions)
				{
					if (regions[j] in good[\"Region\"])
					{
						good[\"Region\"][regions[j]] = 2;
					}
				}
			}
		}
		else
		{
			good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] = 0;
			if (\"detail\" in tableInfo[\"Gene\"][i])
			{
				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
				if (gene != null)
				{
					gene.setAttribute(\"fill-opacity\", \"0.25\");
					gene.setAttribute(\"stroke-opacity\", \"0.25\");
					gene.setAttribute(\"fill\", \"black\");
					gene.setAttribute(\"stroke\", \"none\");
				}
				var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
				if (gene_detail != null)
				{
					gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
					gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
				}
					//Dimming down the genes in preview panes that aren't so annotated
				var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
				if (gene_detail != null)
				{
					gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
					gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
				}
			}
		}
	}

	//Repeating the same for the regions
	for (var i = 0; i < tableInfo[\"Region\"].length; i++)
	{
		if (good[\"Region\"][tableInfo[\"Region\"][i][\"href\"]] != 2)
		{
			good[\"Region\"][tableInfo[\"Region\"][i][\"href\"]] = 0;
			var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
			if (gene != null)
			{
				gene.setAttribute(\"fill-opacity\", \"0.25\");
				gene.setAttribute(\"stroke-opacity\", \"0.25\");
			}

		}
		else
		{
			good[\"Region\"][tableInfo[\"Region\"][i][\"href\"]] = 1;
			var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
			if (gene != null)
			{
				gene.setAttribute(\"stroke\", \"black\");

			}
		}
	}
	
	if (svg_all.hasAttribute(\"tableOn\"))
	{
		makeTable(svg_all.getAttribute(\"tableType\"));
		if (svg_all.hasAttribute(\"expand\"))
		{
			expandTable();
		}
	}
	
}

/*
Function that allows table rows to be selected or deselected
*/
function selectGeneTypeTable(id)
{
	var svg_all = document.getElementById(\"allsvg\");
	var table = null;

	//Getting the table information from the JSON strings
	var tableHead = JSON.parse(tableHeadStr);
	var tableWidth = JSON.parse(tableWidthStr);
	var tableInfo = JSON.parse(tableInfoStr);
	var term2type = JSON.parse(term_2_type);
	
	//Looking to see if a row is already turned on
	if (svg_all.hasAttribute(\"highlight.row.id\"))
	{
		//Get the ID
		var turnOff = document.getElementById(svg_all.getAttribute(\"highlight.row.id\"));
		var turnOffPart = document.getElementById(turnOff.getAttribute(\"href\"));
		if (turnOffPart != null)
		{
			turnSampleOff(turnOffPart);
		}
		svg_all.removeAttribute(\"highlight.row.id\");
		restoreLegend();
		LegendOn=\"yes\";
	}

	//If no row is selected at this point
	if(!(svg_all.hasAttribute(\"Selection\")) || svg_all.getAttribute(\"Selection\") == null)
	{

		//Going through all the genes to check if it has the given id as an annotation
		for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
		{
			//If the gene has type reference, and it contains the search id, flag it by setting good gene as one
			if (\"oth_ref\" in tableInfo[\"Gene\"][i] && tableInfo[\"Gene\"][i][\"oth_ref\"].search(id)>-1)
			{
				good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] = 1;
				if (\"detail\" in tableInfo[\"Gene\"][i])
				{
					//Coloring turned on genes
					var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
					if (gene != null)
					{
						gene.setAttribute(\"fill\", \"blue\");
						gene.setAttribute(\"stroke\", gene.getAttribute(\"fill_old\"));
					}
				}
			}
			else
			{
				//opaquing the genes without the type reference
				if (\"detail\" in tableInfo[\"Gene\"][i])
				{
					var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
					if (gene != null)
					{
						gene.setAttribute(\"fill-opacity\", \"0.25\");
						gene.setAttribute(\"stroke-opacity\", \"0.25\");
						gene.setAttribute(\"fill\", \"black\");
						gene.setAttribute(\"stroke\", \"none\");
					}
					//Dimming down the inner rings if the gene isn't so annotated
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
						gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
					}
					//Dimming down the genes in preview panes that aren't so annotated
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
						gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
					}
				}
			}
		}

		//Repeating the same for the regions
		for (var i = 0; i < tableInfo[\"Region\"].length; i++)
		{
			//Good regions = black outlined
			if (\"oth_ref\" in tableInfo[\"Region\"][i] && (tableInfo[\"Region\"][i][\"oth_ref\"].search(id)>-1))
			{
				var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
				if (gene != null)
				{
					gene.setAttribute(\"stroke\", \"black\");
				}
			}
			else
			{
				//Bad regions = dimmed regions
				var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
				if (gene != null)
				{
					gene.setAttribute(\"fill-opacity\", \"0.25\");
					gene.setAttribute(\"stroke-opacity\", \"0.25\");
				}
			}
		}
		var Leg1 = document.getElementById(\"Legend1\");
		var rects = Leg1.querySelectorAll(\".LegRect\");
		var texts = Leg1.querySelectorAll(\".LegText\");
		var cols = Leg1.querySelectorAll(\".LegCol\");
		for (var i =0; i < rects.length; i++)
		{
			if (term2type[id] == rects[i].getAttribute(\"href\"))
			{
				rects[i].setAttribute(\"opacity\", \"1.0\");
			}
			else
			{
				rects[i].setAttribute(\"opacity\", \"0.25\");
			}
		}
		for (var i= 0; i < texts.length; i++)
		{
			if (term2type[id] == texts[i].getAttribute(\"href\"))
			{
				texts[i].setAttribute(\"opacity\", \"1.0\");
			}
			else
			{
				texts[i].setAttribute(\"opacity\", \"0.25\");
			}
		}
		for (var i =0; i < cols.length; i++)
		{
			if (term2type[id] == cols[i].getAttribute(\"href\"))
			{
				cols[i].setAttribute(\"opacity\", \"1.0\");
			}
			else
			{
				cols[i].setAttribute(\"opacity\", \"0.25\");
			}
		}
		
		

		svg_all.setAttribute(\"Selection\", id);
		LegendOn=\"no\";

	}
	else
	{
		var id_old = svg_all.getAttribute(\"Selection\");
		for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
		{
			//Resetting all the Gene colors + outlines
			if (\"detail\" in tableInfo[\"Gene\"][i])
			{
				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
				if (gene != null)
				{
					gene.setAttribute(\"fill\", gene.getAttribute(\"fill_old\"));
					gene.setAttribute(\"stroke\", \"none\");

					gene.setAttribute(\"fill-opacity\", \"1.0\");
					gene.setAttribute(\"stroke-opacity\", \"1.0\");
				}
				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
				if (gene != null)
				{
					gene.setAttribute(\"fill-opacity\", \"1.0\");
					gene.setAttribute(\"stroke-opacity\", \"1.0\");
				}
				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
				if (gene != null)
				{
					gene.setAttribute(\"fill-opacity\", \"1.0\");
					gene.setAttribute(\"stroke-opacity\", \"1.0\");
				}
			}
		}
		//Resetting all the region colors + outlines

		for (var i = 0; i < tableInfo[\"Region\"].length; i++)
		{
			var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
			if (gene != null)
			{
				gene.setAttribute(\"fill-opacity\", \"1.0\");
				gene.setAttribute(\"stroke\", \"none\");
			}
		}

		//If we are turning things on also (not just off)... repeat the above
		if (id != id_old)
		{
			for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
			{
				if (\"oth_ref\" in tableInfo[\"Gene\"][i] && tableInfo[\"Gene\"][i][\"oth_ref\"].search(id)>-1)
				{
					good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] = 1;
					if (\"detail\" in tableInfo[\"Gene\"][i])
					{
						var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
						if (gene != null)
						{
							gene.setAttribute(\"stroke\", gene.getAttribute(\"fill_old\"));
						}
					}
				}
				else
				{
					if (\"detail\" in tableInfo[\"Gene\"][i])
					{
						var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
						if (gene != null)
						{
							gene.setAttribute(\"fill-opacity\", \"0.25\");
							gene.setAttribute(\"stroke-opacity\", \"0.25\");
							gene.setAttribute(\"fill\", \"black\");
							gene.setAttribute(\"stroke\", \"none\");
						}
						var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
						if (gene_detail != null)
						{
							gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
							gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
						}
						var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
						if (gene_detail != null)
						{
							gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
							gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
						}
					}
				}
			}
			for (var i = 0; i < tableInfo[\"Region\"].length; i++)
			{
				if (\"oth_ref\" in tableInfo[\"Region\"][i] && (tableInfo[\"Region\"][i][\"oth_ref\"].search(id)>-1))
				{
					var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
					if (gene != null)
					{
						gene.setAttribute(\"stroke\", \"black\");
					}
				}
				else
				{
					var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
					if (gene != null)
					{
						gene.setAttribute(\"fill-opacity\", \"0.25\");
						gene.setAttribute(\"stroke-opacity\", \"0.25\");
					}
				}
			}

			svg_all.setAttribute(\"Selection\", id);
			var Leg1 = document.getElementById(\"Legend1\");
			var rects = Leg1.querySelectorAll(\".LegRect\");
			var texts = Leg1.querySelectorAll(\".LegText\");
			var cols = Leg1.querySelectorAll(\".LegCol\");
			for (var i =0; i < rects.length; i++)
			{
				if (term2type[id] == rects[i].getAttribute(\"href\"))
				{
					rects[i].setAttribute(\"opacity\", \"1.0\");
				}
				else
				{
					rects[i].setAttribute(\"opacity\", \"0.25\");
				}
			}
			for (var i= 0; i < texts.length; i++)
			{
				if (term2type[id] == texts[i].getAttribute(\"href\"))
				{
					texts[i].setAttribute(\"opacity\", \"1.0\");
				}
				else
				{
					texts[i].setAttribute(\"opacity\", \"0.25\");
				}
			}
			for (var i =0; i < cols.length; i++)
			{
				if (term2type[id] == cols[i].getAttribute(\"href\"))
				{
					cols[i].setAttribute(\"opacity\", \"1.0\");
				}
				else
				{
					cols[i].setAttribute(\"opacity\", \"0.25\");
				}
			}
			LegendOn=\"no\";

		}
		else
		{
			svg_all.removeAttribute(\"Selection\");
			restoreLegend();
			LegendOn=\"yes\";
		}

	}
}

/*Selects and highlights the regions and genes based on a functional annotation
Highlights genes (inner rings) and regions (outer rings), and trims the tables
If already the function is already selected, it turns it off
*/
function selectGeneType(evt)
{
	var target = evt.target;
	var id = target.getAttribute(\"href\");
	var svg_all = document.getElementById(\"allsvg\");
	var table = null;
	var fillcol = target.getAttribute(\"fillcol\");

	console.log(target);
	//Highlights the function in the legend
	var id_rect = document.getElementById(id + \"rect\");
	if (id_rect != null)
	{
		id_rect.setAttribute(\"visibility\", \"visible\");
	}
	var id_rect2 = document.getElementById(id + \"col\");

	//Gets the table information from the JSON
	var tableHead = JSON.parse(tableHeadStr);
	var tableWidth = JSON.parse(tableWidthStr);
	var tableInfo = JSON.parse(tableInfoStr);
	var term2type = JSON.parse(term_2_type);

	var selectType = \"\";
	var selectCols = \"\";
	var othType = \"\";
	var turnOn = 1;
	
	//Turns off a row highlight if it is on
	if (svg_all.hasAttribute(\"highlight.row.id\"))
	{
		var turnOff = document.getElementById(svg_all.getAttribute(\"highlight.row.id\"));
		turnOff.style.backgroundColor = target.getAttribute(\"old_bg\");
		var turnOffPart = document.getElementById(turnOff.getAttribute(\"href\"));
		if (turnOffPart != null)
		{
			turnSampleOff(turnOffPart);
		}
		svg_all.removeAttribute(\"highlight.row.id\");
	}

	//Highlighing the legend button
	if(!(svg_all.hasAttribute(\"Selection\")))
	{
		svg_all.setAttribute(\"Selection\", \"(\" + id + \")\");
		svg_all.setAttribute(\"SelectCols\", \"(\" + fillcol + \")\");

		id_rect2.style.stroke = \"black\";
		id_rect2.style.strokeWidth = 3;
		selectType =  \"(\" + id + \")\";
		selectCols = \"(\" + fillcol + \")\";

	}
	else
	{
		//Searching for selected type in selection
		selectType = svg_all.getAttribute(\"Selection\");
		var res = selectType.search( \"(\" + id + \")\");
		selectCols = svg_all.getAttribute(\"SelectCols\");
		//If selection has already been selected, turn it off
		if (res >-1)
		{
			selectType = selectType.replace(\"(\" + id + \")\", \"\");
			selectCols = selectCols.replace(\"(\" + fillcol + \")\", \"\");

			id_rect2.style.strokeWidth = 0;
			svg_all.setAttribute(\"Selection\", selectType);
			svg_all.setAttribute(\"SelectCols\", selectCols);

		}
		else
		{
			//else, turn it on (highlight it)

			selectType = selectType +  \"(\" + id + \")\";
			selectCols = selectCols +  \"(\" + fillcol + \")\";

			id_rect2.style.stroke = \"black\";
			id_rect2.style.strokeWidth = 3;
			svg_all.setAttribute(\"Selection\", selectType);
			svg_all.setAttribute(\"SelectCols\", selectCols);
		}
	}

	var ids = selectType.match(/([^\)\(]+)/g); // A list of all the selected ids
	var id_cols = selectCols.match(/([^\)\(]+)/g); // A list of all the selected ids
	
	var kys = Object.keys(good); //array of all the selected genes
	
	if (ids != null)
	{
		for (var i = 0; i < kys.length; i++)
		{
			var kys2 = Object.keys(good[kys[i]]);
			for (var j = 0; j < kys2.length; j++)
			{
				good[kys[i]][kys2[j]] = 0;
			}

		}
		for (var i = 0; i < tableInfo[\"Region\"].length; i++)
		{
			var cnt = 0;
			var genes = reg2gene[tableInfo[\"Region\"][i][\"href\"]].split(\";\");
			var region = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
			
			var fill_val = region.getAttribute(\"fill_old\");
			var addfill = [];
			if (\"type_ref\" in tableInfo[\"Region\"][i])
			{

				for (var j = 0; j < ids.length; j++)
				{
					if (tableInfo[\"Region\"][i][\"type_ref\"].search(ids[j])>-1)
					{
						cnt = cnt + 1;
						addfill.push(id_cols[j]);
						fill_val = id_cols[j];
					}
				}
			}
			if (cnt > 0)
			{
				var findPat = document.getElementById(\"pattern\" + i);
				if (findPat != null)
				{
					findPat.remove();
				}
				var newPat = document.createElementNS(\"http://www.w3.org/2000/svg\", \"pattern\");
				newPat.setAttribute(\"id\", \"pattern\" + i);
				newPat.setAttribute(\"patternUnits\", \"userSpaceOnUse\");
				newPat.setAttribute(\"patternTransform\", \"rotate(\"+ region.getAttribute(\"angle\")+\")\");
				newPat.setAttribute(\"x\", 0);
				newPat.setAttribute(\"y\", 0);
				newPat.setAttribute(\"width\", 10);
				newPat.setAttribute(\"height\", 10);
				newPat.setAttribute(\"viewbox\",\"0 0 10 10\");
				var bar_ht = 10/cnt;
				for (var j = 0; j < cnt; j++)
				{
					var newLine = document.createElementNS(\"http://www.w3.org/2000/svg\", \"rect\");
					newLine.setAttribute(\"x\", 0);
					newLine.setAttribute(\"y\", j * bar_ht);
					newLine.setAttribute(\"width\", 10);
					newLine.setAttribute(\"height\", bar_ht);
					newLine.setAttribute(\"fill\", addfill[j]);
					newPat.appendChild(newLine);
				}
				if (tableInfo[\"Region\"][i][\"Type\"] == \"fGR\")
				{
					for (var j = 0; j < 10; j +=5)
					{
						var newLine = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
						newLine.setAttribute(\"x1\", j);
						newLine.setAttribute(\"y1\", 0);
						newLine.setAttribute(\"x2\", 0);
						newLine.setAttribute(\"y2\", j);
						newLine.setAttribute(\"stroke\", \"black\");
						newLine.setAttribute(\"strokeWidth\", 0.5);
						newPat.appendChild(newLine);
						var newLine2 = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
						newLine2.setAttribute(\"x1\", j);
						newLine2.setAttribute(\"y1\", 10);
						newLine2.setAttribute(\"x2\", 10);
						newLine2.setAttribute(\"y2\", j);
						newLine2.setAttribute(\"stroke\", \"black\");
						newLine2.setAttribute(\"strokeWidth\", 0.5);
						newPat.appendChild(newLine2);
					}
				}
				defs.appendChild(newPat);
				
				fill_val = \"url(#pattern\" + i+\")\";
			}
			region.setAttribute(\"fill\", fill_val);
			region.setAttribute(\"fill-opacity\", 1);

			if (cnt > 0)
			{
				good[\"Region\"][tableInfo[\"Region\"][i][\"href\"]] = 1;
				if (region != null)
				{
				
					region.setAttribute(\"stroke\", \"black\");
					region.setAttribute(\"stroke-width\", \"3\");
					region.setAttribute(\"stroke-opacity\", \"1\");
					
					
					if (tableInfo[\"Region\"][i][\"Type\"] == \"fGR\")
					{
						region.setAttribute(\"fill_prev\", region.getAttribute(\"fill_prev\") + region.getAttribute(\"fill\") + \";\");
						//region.setAttribute(\"stroke-dasharray\", \"1.0,1.0\");
						
					}
					else
					{
						region.setAttribute(\"fill_prev\",  region.getAttribute(\"fill\") + \";\" + region.getAttribute(\"fill_prev\"));
						//region.setAttribute(\"fill-opacity\", 0.15);
						//region.setAttribute(\"stroke-dasharray\", \"1,0\");
					
					}
				}
				for (var j =0; j < genes.length; j++)
				{
					good[\"Gene\"][genes[j]] = 1;
				}

			}
			else
			{
				var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
				good[\"Region\"][tableInfo[\"Region\"][i][\"href\"]] = 0;
				for (var j =0; j < genes.length; j++)
				{
					good[\"Gene\"][genes[j]] = 0;
				}
				if (gene != null)
				{
					region.setAttribute(\"stroke\", null);			
					gene.setAttribute(\"fill-opacity\", \"0.25\");
				}
			}
		}
		for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
		{
			if (\"href\" in tableInfo[\"Gene\"][i])
			{
				if (good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] == 1 && selectType.search(\"(\"+tableInfo[\"Gene\"][i][\"type_ref\"]+\")\")>-1)
				{
					good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] = 1;
					if (\"oth_ref\" in tableInfo[\"Gene\"][i] & othType.search(tableInfo[\"Gene\"][i][\"oth_ref\"]) == -1)
					{
						othType = othType + tableInfo[\"Gene\"][i][\"oth_ref\"];
					}
					var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
					if (gene != null)
					{
						gene.setAttribute(\"fill-opacity\", \"1.0\");
						gene.setAttribute(\"stroke-opacity\", \"1.0\");
						gene.setAttribute(\"stroke\", \"black\");
						gene.setAttribute(\"stroke-width\", \"2.0\");
					}
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\", \"1.0\");
						gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
						gene_detail.setAttribute(\"stroke\", \"black\");
						gene_detail.setAttribute(\"stroke-width\", \"2.0\");
					}
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\", \"1.0\");
						gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
					}
				}

				else
				{
					good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] = 0;
					var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
					if (gene != null)
					{
						gene.setAttribute(\"fill-opacity\", \"0.25\");
						gene.setAttribute(\"stroke-opacity\", \"0.25\");
						gene.setAttribute(\"fill\", \"black\");
						gene.setAttribute(\"stroke\", \"none\");
					}
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
						gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
					}
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
						gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
					}
				}
			}
		}
		for (var i = 0; i < kys.length; i++)
		{
			if (kys[i] != \"Gene\" && kys[i] != \"Region\")
			{
				var kys2 = Object.keys(good[kys[i]]);
				for (var j = 0; j < kys2.length; j++)
				{
					if (kys2[j] in term2type)
					{
						for (j2 = 0; j2 < ids.length; j2++)
						{
							if (term2type[kys2[j]] == ids[j2])
							{
								good[kys[i]][kys2[j]] = 1;
							}
							else
							{
								good[kys[i]][kys2[j]] = 0;
							}
						}
					}
				}
			}
		}
	}
	else
	{
		//Resetting the gene and region colors etc
		var id_old = svg_all.getAttribute(\"Selection\");
		for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
		{
			if (\"detail\" in tableInfo[\"Gene\"][i])
			{
				good[\"Gene\"][tableInfo[\"Gene\"][i][\"href\"]] = 1;

				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
				//Resetting the gene if it exist
				if (gene != null)
				{

					gene.setAttribute(\"fill\", gene.getAttribute(\"fill_old\"));
					gene.setAttribute(\"stroke\", \"none\");

					gene.setAttribute(\"fill-opacity\", \"1.0\");
					gene.setAttribute(\"stroke-opacity\", \"1.0\");
				}
				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
				if (gene != null)
				{
					gene.setAttribute(\"fill-opacity\", \"1.0\");
					gene.setAttribute(\"stroke-opacity\", \"1.0\");
					gene.setAttribute(\"stroke\", \"none\");
				

				}
				var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
				if (gene != null)
				{
					gene.setAttribute(\"fill-opacity\", \"1.0\");
					gene.setAttribute(\"stroke-opacity\", \"1.0\");

				}
			}
		}
		//Resetting all the regions
		for (var i = 0; i < tableInfo[\"Region\"].length; i++)
		{
			var gene = document.getElementById(tableInfo[\"Region\"][i][\"href\"]);
			good[\"Region\"][tableInfo[\"Region\"][i][\"href\"]] = 1;
			if (gene != null)
			{

				gene.setAttribute(\"fill-opacity\", \"1.0\");
				gene.setAttribute(\"fill\", gene.getAttribute(\"fill_old\"));

				gene.removeAttribute(\"stroke\");
			}
		}
		//
		for (var i = 0; i < kys.length; i++)
		{
			if (kys[i] != \"Gene\" && kys[i] != \"Region\")
			{
				var kys2 = Object.keys(good[kys[i]]);
				for (var j = 0; j < kys2.length; j++)
				{
					good[kys[i]][kys2[j]] = 1;
				}
			}
		}

	}
	if (svg_all.hasAttribute(\"tableType\"))
	{
		makeTable(svg_all.getAttribute(\"tableType\"));
	}
}

function expandTable()
{
	var tableDiv = document.getElementById(\"tableDiv\");
	tableDiv.style.top = 0;
	tableDiv.style.left = $border - 10 + \"px\";
	tableDiv.style.height = (2 * ($border + $max_radius)) + \"px\";
	tableDiv.style.width = (2 *($max_radius) + 20) + \"px\";
	tableDiv.style.backgroundColor = \"LightGray\";
	var tableMain = document.getElementById(\"tableMain\");
	tableMain.style.height = (2 * ($border + $max_radius) - 30 )+ \"px\";
	tableMain.style.width = (2 *($max_radius)) + \"px\";
	tableMain.style.left = 10 + \"px\";
	var tableBody = document.getElementById(\"allBody\");
	tableBody.style.height = (2 * ($border + $max_radius) - 60 )+ \"px\";
	tableBody.style.width = (2 *($max_radius)) + \"px\";
	
	var allSVG = document.getElementById(\"allsvg\");
	allSVG.style.opacity = 0.5;
	allSVG.setAttribute(\"expand\", \"yes\");
	var otherTab = document.getElementById(\"otherTable\");
	otherTab.style.width=60 + \"px\"; otherTab.style.height = 30 + \"px\";
	otherTab.style.visibility =\"visible\";
	
	var dropMenu = document.getElementById(\"tableType\");
	dropMenu.style.width = 120 + \"px\";
	dropMenu.style.height = 30 + \"px\";
	dropMenu.style.visibility =\"visible\";
	
	var expandSVG = document.getElementById(\"expandSVG\");
	expandSVG.style.left = (2 *($max_radius)) - 30 + \"px\";
	expandSVG.setAttribute(\"transform\", \"rotate(180)\");
	expandSVG.setAttribute(\"onclick\", \"deflateTable()\");
	
}

function deflateTable()
{
	var tableDiv = document.getElementById(\"tableDiv\");
	tableDiv.style.top = $border + $max_radius - 3.5 * 7.5 + \"px\";  
	tableDiv.style.left = (2.125 * $border) +\"px\";
	tableDiv.style.height = ($max_radius / 4 + 7.5 * 4) + \"px\";
    tableDiv.style.width =  ($max_radius * 2 - 2.25 * $border) + \"px\";
	tableDiv.style.backgroundColor = \"transparent\";

	var tableMain = document.getElementById(\"tableMain\");
	tableMain.style.width = (($max_radius * 2) - (2.25 * $border ))+ \"px\";
	tableMain.style.height = ($max_radius / 4) + \"px\";
	tableMain.style.left = 0 + \"px\";
	var tableBody = document.getElementById(\"allBody\");
	tableBody.style.height = \"$bodyHeight\"; 
	tableBody.style.width = \"100%\";;
	var allSVG = document.getElementById(\"allsvg\");
	allSVG.style.opacity = 1;
	allSVG.removeAttribute(\"expand\");
	var otherTab = document.getElementById(\"otherTable\");
	otherTab.style.width=0 + \"px\"; otherTab.style.height = 0 + \"px\"; otherTab.style.visibility =\"hidden\";
	var dropMenu = document.getElementById(\"tableType\");
	dropMenu.style.width = 0 + \"px\";
	dropMenu.style.height = 0 + \"px\";
	dropMenu.style.visibility =\"hidden\";
	var expandSVG = document.getElementById(\"expandSVG\");
	expandSVG.style.left = (($max_radius * 2) - (2.25 * $border ))- 30 + \"px\";
	expandSVG.setAttribute(\"transform\", \"rotate(0)\");
	expandSVG.setAttribute(\"onclick\", \"expandTable()\");
	
}


/*Change the table from on to off
*/
function changeTableStatus(evt)
{
	var target = document.getElementById(\"ButtonIDTableOn\");
	var svgAll = document.getElementById(\"allsvg\");
	//Turning the table off (ie hiding it and changing the button to \"turn it on\")
	if (svgAll.hasAttribute(\"tableOn\"))
	{
		svgAll.setAttribute(\"tableOff\",svgAll.getAttribute(\"tableOn\"));
		var buttons = document.getElementById(\"tableButtons\");
		buttons.setAttribute(\"visibility\", \"hidden\");
		document.getElementById(\"tableMain\").style.display=\"none\";
		svgAll.removeAttribute(\"tableOn\");
		document.getElementById(\"tableButtonText\").childNodes[0].nodeValue = \"Show Table\";
		target.setAttribute(\"fill\", \"#bbffbb\");
		var expandSVG = document.getElementById(\"expandSVG\");
		expandSVG.style.visibility = \"hidden\";
		
		

	}
	else
	{
		//Turning the table on. If the first time- this is the default

		if (svgAll.hasAttribute(\"tableOff\"))
		{
			svgAll.setAttribute(\"tableOn\",svgAll.getAttribute(\"tableOff\"));
			var buttons = document.getElementById(\"tableButtons\");
			buttons.setAttribute(\"visibility\", \"visible\");
			document.getElementById(\"tableMain\").style.display=\"$disp\";
			svgAll.removeAttribute(\"tableOff\");
			document.getElementById(\"tableButtonText\").childNodes[0].nodeValue = \"Hide Table\";
			target.setAttribute(\"fill\", \"#eeeeff\");
			var expandSVG = document.getElementById(\"expandSVG\");
			expandSVG.style.visibility = \"visible\";
		
		}
		else
		{
			//returning to an already opened table. Just restores

			svgAll.setAttribute(\"tableOn\",\"$default_table\");
			var buttons = document.getElementById(\"tableButtons\");
			buttons.setAttribute(\"visibility\", \"visible\");
			document.getElementById(\"tableMain\").style.display=\"$disp\";
			document.getElementById(\"tableButtonText\").childNodes[0].nodeValue = \"Hide Table\";
			selectTableType(document.getElementById(svgAll.getAttribute(\"tableOn\")), \"$default_table_type\");
			target.setAttribute(\"fill\", \"#eeeeff\");
			var expandSVG = document.getElementById(\"expandSVG\");
			expandSVG.style.visibility = \"visible\";

		}
	}
}

/*Select table wrapper*/
function selectTableButton(evt, id)
{
	var target = evt.target;
	selectTableType(target, id);
}

/* Save the main page image as an SVG or PNG
*/
function saveSVG(fileType)
{
	var event = new Event(\'build\');
	var nw = document.createElement(\"a\");
	var svg_all = document.getElementById(\"allsvg\");
	var coreNum = svg_all.getAttribute(\"coreNum\");
	var svgStr = \"<svg xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" height=\\\"$SVGHEIGHT\\\" width=\\\"$SVGWIDTH\\\">\"+document.getElementById(\"svgOut\"+coreNum).innerHTML +document.getElementById(\"svgArc\"+coreNum).innerHTML+\"</svg>\";
	var inSvg = new Blob([svgStr], {type: \'image/svg+xml\'});
	var svgfile = URL.createObjectURL(inSvg);
	var location = window.location.href;
	var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);

	//Saving the PNG Version
	if (fileType==\"png\")
	{
		var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"svg\");
		newSVG.setAttribute(\"xmlns:xlink\",\"http://www.w3.org/1999/xlink\");
		newSVG.setAttribute(\"height\", \"$SVGHEIGHT\" + \"px\");
		newSVG.setAttribute(\"width\",\"$SVGWIDTH\" +\"px\");
		newSVG.innerHTML = document.getElementById(\"LegendText\").innerHTML + document.getElementById(\"Legend1\").innerHTML+ document.getElementById(\"svgOut\"+coreNum).innerHTML + document.getElementById(\"svgArc\"+coreNum).innerHTML;
		var svgStr = new XMLSerializer().serializeToString(newSVG);

		var canvas = document.createElement(\"canvas\");
		canvas.width = $SVGHEIGHT_BIG;
		canvas.height = $SVGWIDTH_BIG;
		var ctx = canvas.getContext(\"2d\");
		ctx.scale(1, 1);
		
		var img = new Image();
		img.src = \"data:image/svg+xml;base64,\" + window.btoa(svgStr);
		img.onload = function() 
		{
			ctx.drawImage(img, 0, 0, $SVGWIDTH, $SVGHEIGHT, 0, 0, $SVGWIDTH_BIG, $SVGHEIGHT_BIG);
			canvas.style.height = $SVGHEIGHT+\"px\";
			canvas.style.width = $SVGWIDTH+\"px\";
		
			nw.setAttribute(\"download\",\"Main.$filenamePNG\");
			nw.setAttribute(\"href\", canvas.toDataURL(\"image/png\"));
			nw.setAttribute(\"target\", \"_blank\");

			alert(\"Saved PNG to \" + \"$filenamePNG\");
			document.body.appendChild(nw); 
			nw.click();
			document.body.removeChild(nw);
		}
	}
	else	
	{

		nw.setAttribute(\"download\", \"Main.$filenameSVG\");
		nw.setAttribute(\"href\", svgfile);
		nw.setAttribute(\"target\", \"_blank\");
		alert(\"Saved SVG to \" + \"Main.$filenameSVG\");
		document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
	}
}

/*
Save the main page table (assuming that the table is turned on) as a tsv file
*/
function saveTableTxt()
{
	var event = new Event(\'build\');
	var nw = document.createElement(\"a\");
	var txtStr = \"\";
	var tBody = document.getElementById(\"tableMain\").tBodies[0];
	for (var i = 0; i < tBody.rows.length; i++)
	{
		var row = tBody.rows[i];
		txtStr = txtStr + row.cells[0].innerHTML;
		for (var j = 1; j < row.cells.length; j++)
		{
			txtStr = txtStr + \"\t\" + row.cells[j].innerHTML;
		}
		txtStr = txtStr + \"\\n\";
	}
	var inTxt = new Blob([txtStr], {type: \'text\'});
	var txtfile = URL.createObjectURL(inTxt);
	var location = window.location.href;
	var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);

	nw.setAttribute(\"download\", curPath + \"$filenameTXT\");
	nw.setAttribute(\"href\", txtfile);
	nw.setAttribute(\"target\", \"_blank\");
	alert(\"Saved TXT to \" + curPath+ \"$filenameTXT\");
	document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
}

/*
Change the information in the table. Caused by clicking on the buttons surrounding it
*/
function selectTableType(evt, id)
{
	var target = evt;

	//Changing the color of the buttons to show which ones are clicked
	var ButtonID = \"ButtonID\"+id;
	var Button = document.getElementById(ButtonID);
	Button.setAttribute(\"fill\", \"#eeeeff\");
	var arcBorderID = \"ArcButtonID\"+id;
	var arcBorder = document.getElementById(arcBorderID);
	arcBorder.setAttribute(\"stroke\", \"#eeeeff\");

	document.getElementById(\"dropdownID\").innerHTML = id;
	var svgAll = document.getElementById(\"allsvg\");

	//If a table has already been selected, remove previous and set up new one
	if (svgAll.hasAttribute(\"tableType\"))
	{
		var oldButtonID = \"ButtonID\"+svgAll.getAttribute(\"tableType\");
		var oldButton = document.getElementById(oldButtonID);
		var oldArcBorderID = \"ArcButtonID\"+svgAll.getAttribute(\"tableType\");
		var oldArcBorder = document.getElementById(oldArcBorderID);

		//Turning button to off colors/borders
		oldArcBorder.setAttribute(\"stroke\", \"#000000\");
		oldButton.setAttribute(\"fill\", \"#bbffbb\");

		//If a table has already been selected, if the previous and new are not the same
		if (svgAll.getAttribute(\"tableType\") != target.getAttribute(\"id\"))
		{

			//Turning off the row highlight if it is on

			if ((svgAll.hasAttribute(\"highlight.row.id\")))
			{

				var turnOff = document.getElementById(svgAll.getAttribute(\"highlight.row.id\"));
				if (turnOff != null)
				{
					turnOff.style.backgroundColor = target.getAttribute(\"old_bg\");
					var turnOffPart = document.getElementById(turnOff.getAttribute(\"href\"));

					if (turnOffPart != null)
					{
						turnSampleOff(turnOffPart);
					}
					svgAll.removeAttribute(\"highlight.row.id\");
				}
			}


			svgAll.setAttribute(\"tableType\",id);
			var ButtonID = \"ButtonID\"+svgAll.getAttribute(\"tableType\");
			var Button = document.getElementById(ButtonID);

			makeTable(id);
			//Turning button to selected by changing color
			Button.setAttribute(\"fill\", \"#eeeeff\")
		}
	}
	else
	{
		svgAll.setAttribute(\"tableType\",id);
		var ButtonID = \"ButtonID\"+svgAll.getAttribute(\"tableType\");
		var Button = document.getElementById(ButtonID);

		makeTable(id);
		Button.setAttribute(\"fill\", \"#eeeeff\")

	}
}

/* Add the data to the table. Does all the data
*/
function makeTable(id)
{
	var svgAll = document.getElementById(\"allsvg\");
	var table = document.getElementById(\"tableMain\");
	var tableHead = JSON.parse(tableHeadStr);
	var tableWidth = JSON.parse(tableWidthStr);
	var tableInfo = JSON.parse(tableInfoStr);

	var coreNum = null;
	if (svgAll.hasAttribute(\"coreNum\"))
	{
		coreNum = svgAll.hasAttribute(\"coreNum\");
	}
	var select = \"0\";
	if (svgAll.hasAttribute(\"Selection\"))
	{
		select = svgAll.getAttribute(\"Selection\");
	}

	table.innerHTML = \"\";
	var header = document.createElement(\"thead\");
	header.position= \"relative\";
	header.display=\"table\";

	var i; //i is the iterative variable

	//Adding the table header along with the functions that
	for (i =0; i < tableHead[id].length; i++)
	{
		var th = document.createElement(\"th\");
		th.innerHTML = tableHead[id][i];
		th.setAttribute(\"sort_dir\", \"0\");
		th.setAttribute(\"id\", tableHead[id][i]);
		th.setAttribute(\"table_id\",  id+\"Body\");
		var in1 = \"sortColumn(event, \"+i+\")\";
		th.setAttribute(\"onclick\", in1);
		th.setAttribute(\"style\", \"background-color: gray; width: \"+ tableWidth[id][i] + \"%;\");
		header.appendChild(th);
	}

	var body = document.createElement(\"tbody\");
	var cur_clust = svgAll.getAttribute(\"coreNum\");
	var terms = JSON.parse(eval(\"termJSONPlasmid\" + cur_clust));
	//The table style
	body.setAttribute(\"style\", \"overflow-y: scroll; word-wrap: break-word; top: $bodyTop; height: $bodyHeight; width: $bodyWidth; background-color: whitesmoke; display: $disp;\");
	body.setAttribute(\"id\", \"allBody\");
	for (i =0; i < tableInfo[id].length; i++)
	{
			if (good[id][tableInfo[id][i][\"href\"]] == 1)
			{
				var row = document.createElement(\"tr\");
				if (\'href\' in tableInfo[id][i])
				{
					row.setAttribute(\"href\", tableInfo[id][i][\"href\"]);
				}
				if (\'detail\' in tableInfo[id][i])
				{
					row.setAttribute(\"detail\", tableInfo[id][i][\"detail\"]);
			}
				row.setAttribute(\"onclick\", \"runType(\'showRowPlot\', event)\");
				row.setAttribute(\"id\", \"rowID\" + i);
				var j;
				var cellLeft = 0;
				for (j =0; j < tableHead[id].length; j++)
				{
					var td = document.createElement(\"td\");
					var cell = document.createTextNode(tableInfo[id][i][tableHead[id][j]]);
					td.appendChild(cell);
					var newWidth = parseInt(tableWidth[id][j]) * 0.95;
					td.setAttribute(\"style\", \"width: \"+ newWidth.toString() + \"%;display: $disp; left:\" + cellLeft.toString()+\"%;\");
					cellLeft = cellLeft + parseInt(tableWidth[id][j]) * 0.95;
					row.appendChild(td);
				}
				if (\"core_clust\" in tableInfo[id][i] && tableInfo[id][i][\"core_clust\"] == cur_clust)
				{
					body.appendChild(row);
				}
				if (!(\"core_clust\" in tableInfo[id][i]))
				{
					var tmpStr = tableInfo[id][i][\"href\"];
					if (typeof tmpStr != \"undefined\")
					{
						tmpStr = tmpStr.substring(0, tmpStr.length-1);
						if (tmpStr in terms)
						{
							body.appendChild(row);
						}
					}
				}
			}
	}

	//Adding to the table the body and header
	table.appendChild(header);

	table.appendChild(body);



}

/* Wrapper function used show the preview panels from the table
*/
function showPlot(evt) {
	var target = evt.target;
	var svg_all = document.getElementById(\"allsvg\");

	if ((svg_all.hasAttribute(\"highlight.row.id\")))
	{
	
	/*		var turnOff = document.getElementById(svg_all.getAttribute(\"highlight.row.id\"));
			turnOff.style.background = turnOff.getAttribute(\"old_bg\");
			var turnOffPart = document.getElementById(turnOff.getAttribute(\"href\"));
			if (turnOffPart != null)
		{
			turnSampleOff(turnOffPart);
		}
		//svg_all.removeAttribute(\"highlight.row.id\");
	*/
	}
	else
	{
		showSample(target);
	}
}

/* Shows or hides the preview panel in the appropiate corner
*/
function showSample(target)
{
	var id_name = target.getAttribute(\"href\");
	var but_id = target.getAttribute(\"id\");

	//changing the location color to grey to show that it is being shown
	target.setAttribute(\"fill\", \"black\");
	var center = document.getElementById(id_name);

	//Only if there is a panel to show...
	if (center != null)
	{

		var status = center.getAttribute(\"visibility\");
		var svg_all = document.getElementById(\"allsvg\");
		if (svg_all.hasAttribute(\"core_set\"))
		{
			var core_set = svg_all.getAttribute(\"core_set\");
			var core = document.getElementById(core_set);
			turnSampleOff(core);
		}

		//Sets the x and y to the appropiate corner
		var new_x = 0;
		var new_y = 0;
		if (parseFloat(target.getAttribute(\"angle\")) < 90)
		{
			new_x = svg_all.getAttribute(\"width\") - center.getBoundingClientRect().width;
		}
		else if(parseFloat(target.getAttribute(\"angle\")) < 180)
		{
			new_y = svg_all.getAttribute(\"height\") - center.getBoundingClientRect().height;
			new_x = svg_all.getAttribute(\"width\") - center.getBoundingClientRect().width;

		}
		else if(parseFloat(target.getAttribute(\"angle\")) < 270)
		{
			new_y = svg_all.getAttribute(\"height\") - center.getBoundingClientRect().height;

		}

		//Panel is turned off-> need to turn on
		if (status == \"hidden\")
		{
			center.setAttribute(\"visibility\", \"visible\");

			center.setAttribute(\"x\", new_x);
			center.setAttribute(\"y\", new_y);

			svg_all.setAttribute(\"core_set\", id_name);
			svg_all.setAttribute(\"core_button\", but_id);
		}
		else
		{

			//Panel is turned on-> need to turn off
			center.setAttribute(\"visibility\", \"hidden\");
			svg_all.removeAttribute(\"core_set\");
			target.setAttribute(\"fill\", \"#e0e0e0\");
			svg_all.removeAttribute(\"core_button\");
		}
	}
}


/* Closes the preview pane */
function closeCL(evt)
{
	var target = evt.target;
	var id_name = target.getAttribute(\"grp_id\");
	var group = document.getElementById(id_name);
	var svg_all = document.getElementById(\"allsvg\");

	group.setAttribute(\"visibility\", \"hidden\");
	svg_all.removeAttribute(\"core_set\");
	var old_but_id = svg_all.getAttribute(\"core_button\");
	var old_button = document.getElementById(old_but_id);
	var fill_old = old_button.getAttribute(\"fill_old\");
	old_button.setAttribute(\"fill\", fill_old);
	svg_all.removeAttribute(\"core_button\");
}

/*Turns off the menu higlight */
function turnOffHighlight(evt)
{
	var target = evt.target;
	target.style.backgroundColor = \"#bbbbbb\";
}

/* Turns on the the table row highlight */
function turnOnHighlight(evt)
{
	var target = evt.target;
	target.style.backgroundColor = \"#666666\";
}

/*Function that sorts the tables by a column header.
Adds a up or down arrow to the header to show this
*/
function sortColumn(evt, col)
{
	var target = evt.target;
	var table_id = target.getAttribute(\"table_id\");

	//HTML table
	var table = document.getElementById(table_id);

	//Sort direction: ascending or descending
	var dir = target.getAttribute(\"sort_dir\");
	if (table != null)
	{
		var tb = table;
		var tr = Array.prototype.slice.call(tb.rows, 0);
		var i;

		//Is the table already sorted, and by what column
		if (table.hasAttribute(\"SortCol\"))
		{
			var sortCol = table.getAttribute(\"SortCol\");
			var sortColHead = document.getElementById(sortCol);
			sortColHead.innerHTML = sortCol;
			sortColHead.setAttribute(\"sort_dir\", \"0\");
		}

		//If already sorted, change the direction
		if (dir == 1)
		{
			table.setAttribute(\"SortCol\", target.getAttribute(\"id\"));
			target.setAttribute(\"sort_dir\", -1);
			dir = -1;
			target.innerHTML = target.innerHTML.concat(\" \").concat(String.fromCharCode(9650));
		}
		else
		{
			table.setAttribute(\"SortCol\", target.getAttribute(\"id\"));
			target.setAttribute(\"sort_dir\", 1);
			dir = 1;
			target.innerHTML = target.innerHTML.concat(\" \").concat(String.fromCharCode(9660));
		}

		if (isNaN(parseInt(tr[0].cells[col].textContent)))
		{
			tr = tr.sort(function(a, b) { return dir * (a.cells[col].textContent.trim().localeCompare(b.cells[col].textContent.trim()));})
		}
		else
		{
			tr = tr.sort(function(a, b) { return dir * (parseInt(a.cells[col].textContent) - parseInt(b.cells[col].textContent));})
		}

		//Adding the sorted rows to the table body
		for (i = 0; i < tr.length; ++i) {tb.appendChild(tr[i])}
	}
}



/*Shows a preview panel when an appropate row in the table is selected
Also highlights the row (ie turns the background yellow)
*/
function showRowPlot(evt) {
	var target = evt.target;

	//If the target isn't a table row, move up to the parent
	if (target.nodeName.toLowerCase() != \"tr\")
	{
		target = target.parentNode;
	}

	var id_name = target.getAttribute(\"href\");
	var all_table = target.parentNode.parentNode;
	var svg_all = document.getElementById(\"allsvg\");
	var TableType = svg_all.getAttribute(\"tableType\");
	if (target.nodeName.toLowerCase() == \"tr\")
	{
		if (svg_all.hasAttribute(\"core_set\"))
		{
			var core_set = svg_all.getAttribute(\"core_set\");
			var core = document.getElementById(core_set);
			turnSampleOff(core);
		}
	
		//If the row is already highlighted...
		if (target.style.backgroundColor == \"yellow\")
		{
			if ((svg_all.hasAttribute(\"highlight.row.id\")))
			{
				//Turns off the row highlight and closes the panel...
				var turnOff = document.getElementById(svg_all.getAttribute(\"highlight.row.id\"));
				turnOff.style.backgroundColor = target.getAttribute(\"old_bg\");
				var turnOffPart = document.getElementById(target.getAttribute(\"href\"));
				if (turnOffPart != null)
				{
					turnSampleOff(turnOffPart);
				}
			}
			var tableInfo = JSON.parse(tableInfoStr);

			for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
			{
				if (\"detail\" in tableInfo[\"Gene\"][i])
				{
					var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
					if (gene != null)
					{
						gene.setAttribute(\"fill-opacity\", \"1.0\");
						gene.setAttribute(\"stroke-opacity\", \"1.0\");
						gene.setAttribute(\"stroke\", gene.getAttribute(\"fill_old\"));
						gene.setAttribute(\"fill\", gene.getAttribute(\"fill_old\"));
					}
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\", \"1.0\");
						gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
						gene_detail.setAttribute(\"stroke\", gene_detail.getAttribute(\"fill_old\"));
						gene_detail.setAttribute(\"fill\", gene_detail.getAttribute(\"fill_old\"));
					}
					var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
					if (gene_detail != null)
					{
						gene_detail.setAttribute(\"fill-opacity\",\"1.0\");
						gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
					}
					var regionStr = tableInfo[\"Gene\"][i][\"Region\"];
					var regionIds = regionStr.split(\",\");
					for (var j = 0; j < regionIds.length; j++)
					{
						var region = document.getElementById(regionIds[j]);										
						if (region != null && region.hasAttribute(\"stroke\"))
						{
							region.removeAttribute(\"stroke\");
							region.setAttribute(\"fill\", region.getAttribute(\"fill_old\"));
						}
					}
				}
			}
					
			
			if (TableType != \"Gene\" && TableType != \"Region\")
			{
				selectGeneTypeTable(id_name);

			}
			
			
			svg_all.removeAttribute(\"muted_genes\")
		
			svg_all.removeAttribute(\"highlight.row.id\")
		}
		else
		{
			//Sets the old_background variable so it can be restored in time...
			target.setAttribute(\"old_bg\", target.style.backgroundColor);
			target.style.backgroundColor = \"yellow\";

			//Checking to see if row is already on. If so turns it \"off\"
			if ((svg_all.hasAttribute(\"highlight.row.id\")))
			{
				var turnOff = document.getElementById(svg_all.getAttribute(\"highlight.row.id\"));
				turnOff.style.backgroundColor = turnOff.getAttribute(\"old_bg\");
				var turnOffPart = document.getElementById(turnOff.getAttribute(\"href\"));
				if (turnOffPart != null)
				{
					turnSampleOff(turnOffPart);
				}
				if (svg_all.hasAttribute(\"old_gene\"))
				{
					var old_gene = document.getElementById(svg_all.getAttribute(\"old_gene\"));
					old_gene.setAttribute(\"stroke\", \"black\");
					old_gene.setAttribute(\"stroke-width\", \"1\");
					var detail = svg_all.getAttribute(\"old_gene\") + \"detail\";
					var old_gene = document.getElementById(detail);
					old_gene.setAttribute(\"stroke\", \"black\");
					old_gene.setAttribute(\"stroke-width\", \"1\");

				}
			}
			if (svg_all.hasAttribute(\"muted_genes\"))
			{
				var tableInfo = JSON.parse(tableInfoStr);
				for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
				{
					if (\"detail\" in tableInfo[\"Gene\"][i])
					{
						var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
						if (gene != null)
						{
							gene.setAttribute(\"fill-opacity\", \"1.0\");
							gene.setAttribute(\"stroke-opacity\", \"1.0\");
							gene.setAttribute(\"stroke\", gene.getAttribute(\"fill_old\"));
							gene.setAttribute(\"fill\", gene.getAttribute(\"fill_old\"));
						}
						var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
						if (gene_detail != null)
						{
							gene_detail.setAttribute(\"fill-opacity\", \"1.0\");
							gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
							gene_detail.setAttribute(\"stroke\",\"none\");
							gene_detail.setAttribute(\"fill\", gene_detail.getAttribute(\"fill_old\"));
						}
						var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
						if (gene_detail != null)
						{
							gene_detail.setAttribute(\"fill-opacity\", \"1.0\");
							gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
						}
						var region = document.getElementById(tableInfo[\"Gene\"][i][\"Region\"]);
						if (region != null && region.hasAttribute(\"stroke\"))
						{
							region.removeAttribute(\"stroke\");
							region.setAttribute(\"fill\", region.getAttribute(\"fill_old\"));
						}
					}
				}
					
			}
			var id = document.getElementById(id_name);
			if (TableType == \"Gene\" || TableType == \"Region\")
			{
				if (TableType == \"Gene\")
				{
					if (!svg_all.hasAttribute(\"muted_genes\") || svg_all.getAttribute(\"muted_genes\") != id_name)
					{
						var tableInfo = JSON.parse(tableInfoStr);
						svg_all.setAttribute(\"muted_genes\", id_name);
						var hits = 0;
						for (var i = 0; i < tableInfo[\"Gene\"].length; i++)
						{
							if (\"href\" in tableInfo[\"Gene\"][i])
							{
								if (tableInfo[\"Gene\"][i][\"href\"] == id_name)
								{
									hits = 1;
									var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
									if (gene != null)
									{
										gene.setAttribute(\"fill-opacity\", \"1.0\");
										gene.setAttribute(\"stroke-opacity\", \"1.0\");
										gene.setAttribute(\"stroke\", gene.getAttribute(\"fill_old\"));
										gene.setAttribute(\"fill\", gene.getAttribute(\"fill_old\"));
									}
									var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
									if (gene_detail != null)
									{
										gene_detail.setAttribute(\"fill-opacity\", \"1.0\");
										gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
										gene_detail.setAttribute(\"stroke\", \"black\");
										gene_detail.setAttribute(\"stroke-width\", \"2.0\");
									}
									var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
									if (gene_detail != null)
									{
										gene_detail.setAttribute(\"fill-opacity\", \"1.0\");
										gene_detail.setAttribute(\"stroke-opacity\", \"1.0\");
										

									}
									var regionStr = tableInfo[\"Gene\"][i][\"Region\"];
									var regionIds = regionStr.split(\",\");
									for (var j = 0; j < regionIds.length; j++)
									{
										var region = document.getElementById(regionIds[j]);
									
										if (region != null)
										{
											region.setAttribute(\"stroke\", \"black\");;
											region.setAttribute(\"fill\", \"gray\");
										}
									}
								}
								else
								{
									var gene = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]);
									if (gene != null)
									{
										gene.setAttribute(\"fill-opacity\", \"0.25\");
										gene.setAttribute(\"stroke-opacity\", \"0.25\");
										gene.setAttribute(\"fill\", \"black\");
										gene.setAttribute(\"stroke\", \"none\");
									}
									var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detail\");
									if (gene_detail != null)
									{
										gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
										gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
									}
									var gene_detail = document.getElementById(tableInfo[\"Gene\"][i][\"detail\"]+\"detailHT\");
									if (gene_detail != null)
									{
										gene_detail.setAttribute(\"fill-opacity\", \"0.25\");
										gene_detail.setAttribute(\"stroke-opacity\", \"0.25\");
									}
								}
						
							}
						}
						if (hits == 0)
						{
							var region = document.getElementById(id_name);
							if (region != null)
							{
								region.setAttribute(\"stroke\", \"black\");;
								region.setAttribute(\"fill\", \"gray\");
							}		
							
						}
					}
					else
					{
						svg_all.removeAttribute(\"muted_genes\");
					}
				}
				if (TableType == \"Region\")
				{					
					showSample(id);
				}
				if (target.hasAttribute(\"detail\"))
				{
					var detail = target.getAttribute(\"detail\");
					var gene = document.getElementById(detail);
					if (gene != null)
					{
						gene.setAttribute(\"stroke\", \"Yellow\");
						gene.setAttribute(\"stroke-width\", \"3\");
						detail = detail + \"detail\";
						var gene = document.getElementById(detail);
						gene.setAttribute(\"stroke\", \"Red\");
						gene.setAttribute(\"stroke-width\", \"3\");

						svg_all.setAttribute(\"old_gene\", target.getAttribute(\"detail\"));
					}
				}
			}
			else
			{
				selectGeneTypeTable(id_name);
			}
			svg_all.setAttribute(\"highlight.row.id\", target.getAttribute(\"id\"));
		}
	}
}

/*Makes the gene page for plasmid sequences. */
function writeFasta(id,plasmid)
{
	var fasta = JSON.parse(eval(\"fastaJSON\"+plasmid));
	var type = JSON.parse(outJSON);
	var terms = JSON.parse(eval(\"termJSON\"+plasmid));
		//If the gene passed along ID, write the html
	if (id in fasta)
	{
		var html = \"<html><script type=\'text/javascript\' src=\'./json/\"+id+\".allFasta.json\'></script><script type=\'text/javascript\' src=\'https\://s3.eu-central-1.amazonaws.com/cdn.bio.sh/msa/latest/msa.min.gz.js\'></script>\";
		if (treeJSON != null)
		{
			html = html + \"<script type=\'text/javascript\' src=\'./json/tree.json\'></script>\";
		}
				html = html + \"<script type=\'text/javascript\' src=\'./scripts/geneFile.functions.js\'></script><body style=\'font-family:Courier New\'>\";
				//Adding tab radio buttons
				html = html + \"<div id=\'select\'><form><input type=\'radio\' name=\'display\' value=\'summary\' checked onclick=\'changeDiv(event)\'>Summary</input>\";
				html = html + \"<input type=\'radio\' name=\'display\' value=\'centFasta\' onclick=\'changeDiv(event)\'>Centroid Fasta</input>\";
				html = html + \"<input type=\'radio\' name=\'display\' value=\'multFastaView\' onclick=\'changeDiv(event)\'>Multi-Fasta Viewer</input>\";
				html = html + \"<input type=\'radio\' name=\'display\' value=\'allFasta\' onclick=\'changeDiv(event)\'>Download Multi-Fasta File</input>\";
				if (treeJSON != null)
				{
					html = html + \"<input type=\'radio\' name=\'display\' value=\'tree\' onclick=\'changeDiv(event)\'>View Phylogeny</input>\";

				}
				//Adding in the summary information in a table
				html = html + \"</form></div><div id=\'summary\' style=\'visibility:visible;\'><table >\";
				for (var j in type)
				{
					var i = type[j];
					html = html + \"<tr><td style=\'border:1px solid black;\' id=\'ID\" + i + \"\'>\"+i+\"</td>\";
					var outVar = fasta[id][i];
					if (i == \"Sequence\")
					{
						outVar = outVar.replace(/--ret--/g, \"<br>\");
					}
					if (i == \"Associated Terms\")
					{
						var tmpVar = outVar.split(\";\");
						outVar = \"\";
						for (j in tmpVar)
						{
							if (tmpVar[j] in terms)
							{
								outVar = outVar + tmpVar[j] +\"(\"+ terms[tmpVar[j]] +\")<br>\";
							}
						}
					}
					html = html + \"<td style=\'border:1px solid black;\' id=\'\" + i + \"Value\'>\"+outVar+\"</td></tr>\";
				}
				html = html + \"</table></div><div id=\'centFasta\' style=\'visibility:hidden;\'>\";
				//adding in the fasta pages
				var outVar = fasta[id][\"Sequence\"];
				outVar = outVar.replace(/--ret--/g, \"<br>\");
				var fastaHeader = \"&gt;\";
				fastaHeader = fastaHeader + \"Centroid_\" + id + \"<br>\";
				html = html + fastaHeader + outVar + \"</p></div><div id=\'multFastaView\' style=\'visibility:hidden;\'><div id=\'multFastaViewInner\'></div>Viewer is supported by <a href=\'msa.biojs.net\'>MSA Viewer</a></div><div id=\'allFasta\' style=\'visibility:hidden;\'>\";
				//Making tree page if nessicary
				if (treeJSON != null)
				{
					html = html + \"</div><div id='tree' style=\'visibility:hidden;\'>\";

					html = html + \"<table ><tr height=\\\"\" + (tree_size * 0.2) + \"\\\"><td style=\\\"vertical-align:top;\\\" rowspan=\\\"3\\\">\"

					html = html + \"<div><svg version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" id=\\\"tree_svg\\\" height=\\\"\"+ tree_size+\"\\\" width=\\\"\"+ tree_size+\"\\\"></svg>\";

					html = html + \"</div></td><td style=\\\"vertical-align:top;\\\"><div><h3>Tree Functions</h3>Outer Ring Tree Level:<button onclick=\\\"changeTree(\'-\')\\\">-</button><b id = \\\"LevelCount\\\">4</b><button onclick=\\\"changeTree(\'+\')\\\">+</button>\";

					html = html + \"<form>Phylogeny Style:<input type=\\\"radio\\\" name=\\\"type\\\" value=\\\"Circular\\\" checked=\\\"true\\\" onclick=\\\"changeType(event)\\\"/>Circular<input type=\\\"radio\\\" name=\\\"type\\\" value=\\\"Linear\\\" onclick=\\\"changeType(event)\\\"/>Linear</form>\";

					html = html + \"<svg version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" id=\\\"save_svg\\\" height=\\\"\"+ (tree_size *0.1)+\"\\\" width=\\\"\"+ (tree_size*0.4)+\"\\\"><text x=\\\"2\\\" y=\\\"14\\\">Node Colors:</text><circle cx=\\\"155\\\" r=\\\"4\\\" cy=\\\"8\\\" fill=\\\"black\\\"/><text x=\\\"162\\\" y=\\\"14\\\">Present</text><circle cx=\\\"245\\\" r=\\\"4\\\" cy=\\\"8\\\" fill=\\\"LightGray\\\"/><text x=\\\"252\\\" y=\\\"14\\\">Missing</text>$disk_image_svg $disk_image_png</svg></div></td></tr><tr height=\\\"\"+ (tree_size *0.4)+\"\\\" >\";
					html = html + \"<td style=\\\"vertical-align:top;\\\"><div id=\\\"gene_list\\\" style=\\\"height:\"+ (tree_size*0.4)+\"; overflow-y: auto;\\\"></div></td></tr><tr height=\\\"\"+ (tree_size*0.15)+\"\\\"><td style=\\\"vertical-align:bottom;\\\">\";
					html = html + \"<div id=\\\"percentage_image\\\"></div></td></tr><tr height=\\\"\"+ (tree_size*0.25)+\"\\\">\";
					html = html + \"<td><table><tr><td id=\\\"SelectMetatype\\\" width=\\\"\"+ (tree_size*0.4)+\"\\\" style=\\\"vertical-align:top;\\\"><h3>Metadata Type</h3>\";
					for (var a in metadata_head)
					{
						html = html + \"<input type=\\\"checkbox\\\" id=\\\"select\"+a+\"\\\" value=\\\"\" + a+ \"\\\" onchange=\\\"change_metadata_status(event)\\\" />\"+a+\"<p>\";
					}
					html = html +\"</td><td id = \\\"Legend\\\" style=\\\"vertical-align:top;\\\"><h3>Legend</h3><svg id=\\\"LegendSVG\\\" version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" height=\\\"\"+ (tree_size*0.2)+\"\\\" width=\\\"\"+ (tree_size*0.6)+\"\\\"></svg></td></tr></table><td style=\\\"vertical-align:top;\\\"><h3>Percentages</h3><div id =\\\"percDiv\\\"></div></td></tr></table>\";
				}
				html = html +\"</div></body></html>\";
				//Writing html to a new page
				var newWind = window.open();
				newWind.document.write(html);
				newWind.document.close();

			}
			else
			{
				alert(\"No Sequence information for \" + id);
			}
}

var outJSON=\'[\"Cluster\", \"Name\", \"Region\", \"# of Genomes\",\"Functional IDs\",\"Associated Terms\", \"Mean\", \"Standard Deviation\", \"Minimum Length\", \"Maximum Length\", \"Sequence\"]\';

";

}

#Getting the ontology from the function files info
sub load_functions {

    my $funcStr = $_[0];
    my ( $funcID, $input_file, $obo_file, $names, $cols ) = @$funcStr
      or die("Incorrect Function Input for $funcStr. Please check pan_html.ph file or manual");

    #Opening the GO
    open( GO, "<", $input_file )
      or die("Cannot not find $input_file from Configure file. Please check $func_file\n\n")
      ;    #mapterm file from annotation pipeline
    my $isNameFile = 0;

    #Checking if the names file exist
    if ( -e $names ) {

        $isNameFile = 1;

    }

    my $id;

    #Checking for ontology files
    if ( -e $obo_file ) {

        warn "Reading in OBO ontology file $obo_file...\n";
        open( OBO, "<", $obo_file ) or die("Cannot open obo file");
        while ( !eof(OBO) ) {

            my $v = <OBO>;
            chomp($v);
            if ( $v =~ /\A([^:]+):(\s+)([^\n]+)/ ) {

                my $type = $1;
                my $val  = $3;
                $val =~ s/\;//g;
                $val =~ s/\'//g;
                $val =~ s/\"//g;
                if ( $type eq "id" ) { $id = $val; }
                if ( $type eq "name" && $id ) { $ont{$id} = $val; $id = ""; }

            }

        }
        $head{$funcID}      = [ "ID", "Name", "# of Genes with Hits" ];
        $col_width{$funcID} = [ .20,  .65,    .15 ];

    }

    #if the color is a file, open it and read it
    if ( -e $cols ) {

        open( COLS, "<", $cols );
        while ( !eof(COLS) ) {

            my $in = <COLS>;
            chomp($in);
            my @in = split "\t", $in;
            $rgb{ $in[0] } = $in[1];

        }

    } else {

        #Otherwise, treat it as an rgb string
        $rgb{$names} = $cols;

    }

    #Reading in the function map file
    warn "Reading in Function file $input_file ...\n";
    my %gene2;
    my $in = <GO>;
    chomp($in);
    my @head = split "\t", $in;

    #This is the header
    for ( my $i = 1 ; $i < scalar(@head) ; $i++ ) {

        $head[$i] =~ s/\s/\_/g;

    }

    while ( !eof(GO) ) {

        my $in = <GO>;
        chomp($in);

        #cleaning up the string
        $in =~ s/\;//g;
        $in =~ s/\'//g;
        $in =~ s/\"//g;

        my @in_list = split "\t", $in;

        my $id = $in_list[0];
        my $c  = 0;
        if ( $ont{$funcID}->{$id} ) { $c = scalar( @{ $ont{$funcID}->{$id} } ) }

        #Gets the id/gene number for the gene
        $id =~ /centroid_(\d+)/;
        my $id_n = $1;

        #Making sure that the gene exists before assigning the functional information
        if ( $gene_num{$id_n} ) {

            #going through the different ids
            while ( $in_list[1] =~ /([A-Z0-9\:]+)/g ) {

                my $go_id = $1;    #This is the sub-function id
                if ( $go_name{$funcID}->{$go_id} && $go_name{$funcID}->{$go_id} !~ /"CL_$id_n;/ ) {

                    #
                    $jso{$funcID}->[ $go_name{$funcID}->{$go_id} - 1 ]->{ $head{$funcID}->[2] }++;
                    $gene2{$id_n}->{$funcID}->{$go_id} = 1;
					
                } else {

                    if ( %go_name && $go_name{$funcID} ) {

                        $go_name{$funcID}->{$go_id} = 1 + scalar( keys( %{ $go_name{$funcID} } ) );

                    } else {

                        $go_name{$funcID}->{$go_id} = 1;

                    }
                    if ( $head{$funcID} ) {

                        $jso{$funcID}->[ $go_name{$funcID}->{$go_id} - 1 ]->{"href"}                = "$go_id;";
                        $jso{$funcID}->[ $go_name{$funcID}->{$go_id} - 1 ]->{ $head{$funcID}->[0] } = "$go_id";
                        $jso{$funcID}->[ $go_name{$funcID}->{$go_id} - 1 ]->{ $head{$funcID}->[1] } = $ont{$go_id};
                        $jso{$funcID}->[ $go_name{$funcID}->{$go_id} - 1 ]->{ $head{$funcID}->[2] }++;

                    }
                    $gene2{$id_n}->{$funcID}->{$go_id} = 1;

                }

                my $tmp_go = "$go_id;";

                #Checks if the sub-function it is in the list or the if the sub-function exists
                if ( !$jso{Gene}->[ $gene_num{$id_n} ]->{oth_ref} ) {

                    $jso{Gene}->[ $gene_num{$id_n} ]->{oth_ref} = $tmp_go;

                } elsif ( $jso{Gene}->[ $gene_num{$id_n} ]->{oth_ref} !~ /$tmp_go/ ) {

                    $jso{Gene}->[ $gene_num{$id_n} ]->{oth_ref} .= $tmp_go;

                }

                #If not in a name file, just assign the name in the string to all genes
                if ( !$isNameFile ) {

                    if ( $gene_num{$id_n} && $jso{Gene}->[ $gene_num{$id_n} ]->{type_ref} eq $defaultFunct ) {

                        $gene_func{ "CL_" . $id_n } = $names;
                        $func_list{$names} = 1;
                        if ( $head{$funcID} ) {

                            $jso{$funcID}->[ $go_name{$funcID}->{$go_id} - 1 ]->{type_ref} = $names;

                        }
                        $jso{Gene}->[ $gene_num{$id_n} ]->{type_ref} = $names;
                        $jso{Gene}->[ $gene_num{$id_n} ]->{ $head{Gene}->[4] } = $names;

                    }

                }

            }

        }

    }

    close(GO);

    #Getting the type reference name (used in legend)
    if ($isNameFile) {

        open( FUNC, "<", $names ) or die();    #map of go term to cluster ID

        while ( !eof(FUNC) ) {

            my $v = <FUNC>;
            chomp($v);
            my @n = split "\t", $v;
            $n[0] =~ /_(\d+)/;
            my $num = $1;
            if ( $gene_num{$num} && $jso{Gene}->[ $gene_num{$num} ]->{type_ref} eq $defaultFunct ) {

                $jso{Gene}->[ $gene_num{$num} ]->{type_ref} = $n[1];
                $jso{Gene}->[ $gene_num{$num} ]->{ $head{Gene}->[4] } = $n[1];
                if ( scalar(@n) == 3 && $head{$funcID} && $go_name{$funcID}->{ $n[2] } ) {

                    $jso{$funcID}->[ $go_name{$funcID}->{ $n[2] } - 1 ]->{"type_ref"} = $n[1];

                }
                $func_list{ $n[1] } = 1;
                $gene_func{ "CL_" . $num } = $n[1];

            }

        }

    }
	

}

#Reading in the fasta files and the PanACEA Format Files
sub get_gene_info_new {

    while ( $mult_align_dir =~ /([^,]+)/g ) {

        my $try_dir = $1;

        #Looking for multiple alignment directory...
        if ( -e $try_dir && -d $try_dir ) {

            warn "Found multiple sequence alignment directory. Using these as the default seqeunces....\n\n\n";
            opendir( AFA_DIR, $try_dir );
            my @files = readdir(AFA_DIR);

            #Probably allow for multiple fasta names
            foreach my $afa_file (@files) {

                if ( $afa_file =~ /(\S+).(afa|fasta|fa|mfa|nuc)\Z/ ) {

                    #Reading in sequences into seq
                    open( MSA_FASTA, "<", "$try_dir/$afa_file" );
                    while ( !eof(MSA_FASTA) ) {

                        my $in = <MSA_FASTA>;
                        chomp($in);
                        while ( $in =~ />(\S+)/ && !eof(MSA_FASTA) ) {

                            my $id = $1;
                            $in = <MSA_FASTA>;
                            my $tmpSeq;

                            #Adding --ret-- to end of the sequence
                            while ( !eof(MSA_FASTA) && $in !~ />/ ) {

                                chomp($in);
                                $tmpSeq .= $in . "--ret--";
                                $in = <MSA_FASTA>;

                            }
                            if ( eof(MSA_FASTA) ) {

                                chomp($in);
                                $tmpSeq .= $in . "--ret--";

                            }

                            #adding of sequence to $seqs to seq and the length to len
                            if ( !$seqs{all}->{$id}->{seq} ) {

                                $seqs{all}->{$id}->{seq} = $tmpSeq;
                                $tmpSeq =~ s/--ret--//g;
                                $tmpSeq =~ s/-//g;
                                $seqs{all}->{$id}->{len} = length($tmpSeq);

                            }

                        }

                    }

                }

            }

        }

        #Looking for the multiple alignment file (not directory)
        if ( -e $try_dir && -f $try_dir ) {

            open( MSA_FASTA, "<", "$try_dir" );
            while ( !eof(MSA_FASTA) ) {

                #
                my $in = <MSA_FASTA>;
                chomp($in);
                while ( $in =~ />(\S+)/ && !eof(MSA_FASTA) ) {

                    my $id = $1;
                    $in = <MSA_FASTA>;
                    my $tmpSeq;
                    while ( !eof(MSA_FASTA) && $in !~ />/ ) {

                        chomp($in);
                        $tmpSeq .= $in . "--ret--";
                        $in = <MSA_FASTA>;

                    }
                    if ( eof(MSA_FASTA) ) {

                        chomp($in);
                        $tmpSeq .= $in . "--ret--";

                    }
                    if ( !$seqs{all}->{$id}->{seq} ) {

                        $seqs{all}->{$id}->{seq} = $tmpSeq;
                        $tmpSeq =~ s/--ret--//g;
                        $tmpSeq =~ s/-//g;
                        $seqs{all}->{$id}->{len} = length($tmpSeq);

                    }

                }

            }

        }

    }
    warn "Reading in the PanACRA Flat File Info...\n";

    #What it says above...
    open( PANACEA_FF, "<", $file_in1 ) or die("Cannot read in PanACEA Flat File. Please see manual for more information");
    my $c = 0;
    while ( !eof(PANACEA_FF) ) {

        my $in = <PANACEA_FF>;

        #PanACEA Flat Files works modularly, with each module starting with START following by TYPE then the ID
        while ( $in && $in =~ /\ASTART\t(\S+)\t(\S+)/ && !eof(PANACEA_FF) ) {

            my $type = $1;
            my $id   = $2;
            my $tmp;    #stores the information for each module in a hash
            $in = <PANACEA_FF>;    #Reading in the next line (to start the loop)
            $in =~ /\A([^\n\r]+)/;
            $in = $1;

            #Looping to and adding in
            while ( $in !~ /\AEND\Z/ && !eof(PANACEA_FF) ) {

                my @n = split "\t", $in; #each line has the variable id in column 1 and the value in column 2 (tab seperated)
                $tmp->{ $n[0] } = $n[1];
                $in = <PANACEA_FF>;
                $in =~ /\A([^\n\r]+)/;
                $in = $1;

            }

            #After reaching END, the module is finished. Now time to store the variables in the appropiate variable

            #Reading in chrom types
            if ( $type eq "CHROM" && $tmp->{is_core} == 1 ) {

                $cores{$id}->{type} = $tmp->{type};
                $cores{$id}->{sz}   = $tmp->{sz};

            }

            #If it is core gene, push variables into $core_list and $core_size
            if ( $type eq "CORE" ) {

                $tmp->{len} = 1;
                my @in = ( $tmp->{chr}, $id, $tmp->{start}, $tmp->{end}, $tmp->{def}, $tmp->{type}, $tmp->{len} );
                push @{ $core_list{ $tmp->{chr} } }, \@in;
                if ( $tmp->{type} eq "CL" ) {

                    $core_size{$id} = $tmp->{end} - $tmp->{start};

                }

            }

#If it is fgr region (ie a psuedo-core gene), push variables into $list. Also gets the list of genes (gen), the bounding genes and the genomes
            if ( $type eq "FGR" ) {

                $list{$id}->{ $tmp->{order} }->{st} = "";
                if ( $tmp->{st} ) {

                    $list{$id}->{ $tmp->{order} }->{st} = $tmp->{st};    #starting bounding gene

                }
                $list{$id}->{ $tmp->{order} }->{end} = "";
                if ( $tmp->{end} ) {

                    $list{$id}->{ $tmp->{order} }->{end} = $tmp->{end};    #ending bounding gene

                }
                $list{$id}->{ $tmp->{order} }->{gen} = $tmp->{gen};                 #list of genes
                $list{$id}->{ $tmp->{order} }->{cnt} = () = $tmp->{gen} =~ /:/g;    #Number of genes
                my @arr = split ":", $tmp->{order};
                if ( $tmp->{end} ) {

                    pop @arr;

                }
                if ( $tmp->{st} ) {

                    shift @arr;

                }
                $list{$id}->{ $tmp->{order} }->{arr} = \@arr;                       #Array of genes in order
                foreach my $gene (@arr) {

                    $fgi_member{$gene} = $id;    #Adding a map of the gene name to the fGR id

                }

            }

            #Non-core chromosomes, ie circular fGR ie plasmid-like structures
            if ( $type eq "CHROM" && $tmp->{is_core} == 0 ) {

                $fgi_grps{$id}->{type} = $tmp->{type};
                $fgi_grps{$id}->{sz}   = $tmp->{sz};

            }

            #Getting cluster/gene information
            if ( $type eq "CLUSTER" ) {

                $clusters{$id}->{protein_name}   = $tmp->{protein_name};      #Name
                $clusters{$id}->{num_of_members} = $tmp->{num_of_members};    #Number of genomes that contain the cluster
                $clusters{$id}->{centroid}       = $tmp->{centroid};          #ID of the centroid member
                if ( $tmp->{start} ) {

                    $clusters{$id}->{start} = $tmp->{start};                  #BP of gene start
                    $clusters{$id}->{end}   = $tmp->{end};                    #BP of gene end

                } else {

                    $clusters{$id}->{start} = 1;
                    $clusters{$id}->{end}   = 1;

                }

                #Getting the genomes & sequences for the genes
                #Also getting the length statistics: min, max, mean and standard deviation
                my @genomes = split ";", $tmp->{genomes};
                my @names   = split ";", $tmp->{names};
                my $minLen  = -1;
                my $maxLen   = -1;    #default max and min lengths
                my $sum_len  = 0;
                my $sum_len2 = 0;     #sum of lengths and sum of square of sequence lengths
                my $n        = 0;     #number of sequences
                for ( my $i = 0 ; $i < scalar(@genomes) ; $i++ ) {

                    if ( $seqs{all}->{ $names[$i] } ) {

                        my $len = $seqs{all}->{ $names[$i] }->{len};
                        if ( $seqs{all}->{ $names[$i] }->{seq} ) {

                            $seqs{byClust}->{$id}->{ $genomes[$i] } = $seqs{all}->{ $names[$i] }->{seq};

                        } else {

                            $seqs{byClust}->{$id}->{ $genomes[$i] } = "NA";

                        }
                        $seqs{byClustId}->{$id}->{ $genomes[$i] } = $names[$i];
                        if ( $minLen == -1 || $minLen > $len ) {

                            $minLen = $len;

                        }
                        if ( $maxLen == -1 || $maxLen < $len ) {

                            $maxLen = $len;

                        }
                        $sum_len  += $len;
                        $sum_len2 += $len * $len;
                        $n++;

                    } else {

                        warn "No sequence for ", $names[$i], "\n";

                    }

                }
                $c++;
                if ( $jso{Gene} ) { $c = scalar( @{ $jso{Gene} } ); }    #$c is the count number for the jso array gene
                $gene_num{$id} = $c;

                #Calculating the mean (as an int) and standard deviation for each
                #stores both in the gene info length and the appropiate table JSON row
                if ( $n > 0 ) {

                    $jso{Gene}->[$c]->{Mean} = int( 0.5 + $sum_len / $n );
                    $geneLenInfo{$id}->{mean} = int( 0.5 + $sum_len / $n );
                    $jso{Gene}->[$c]->{"Standard Deviation"} =
                      int( 0.5 + 100 * sqrt( $sum_len2 / $n - ( $sum_len / $n ) * ( $sum_len / $n ) ) ) / 100;
                    $geneLenInfo{$id}->{sd} =
                      int( 0.5 + 100 * sqrt( $sum_len2 / $n - ( $sum_len / $n ) * ( $sum_len / $n ) ) ) / 100;
                    $jso{Gene}->[$c]->{"Maximum Length"} = $maxLen;
                    $geneLenInfo{$id}->{maxLen}          = $maxLen;
                    $jso{Gene}->[$c]->{"Minimum Length"} = $minLen;
                    $geneLenInfo{$id}->{minLen}          = $minLen;

                } else {

                    warn "Could not find any sequences for $id...\n";

                }

                #Adding the appropiate gene info to the appropiate column for the table for the gene row
                $jso{Gene}->[$c]->{ $head{Gene}->[0] } = $id;
                $jso{Gene}->[$c]->{ $head{Gene}->[1] } = $tmp->{protein_name};
                $jso{Gene}->[$c]->{ $head{Gene}->[2] } = $tmp->{num_of_members};
                $jso{Gene}->[$c]->{ $head{Gene}->[3] } = "";
                $jso{Gene}->[$c]->{ $head{Gene}->[4] } = $defaultFunct;

                #more table info that is used to choose rows
                $jso{Gene}->[$c]->{type_ref} = $defaultFunct;
                $jso{Gene}->[$c]->{oth_ref}  = "";
                $gene_func{$id}              = $defaultFunct;

                #getting the centroid sequence
                $seqs{cent}->{$id}->{seq} = $seqs{all}->{ $clusters{$id}->{centroid} }->{seq};

                #Tries to get the number of genomes in the total. Probably should be the scalar(@head)-7;
                if ( $GENOME_NUM < $clusters{$id}->{num_of_members} ) {

                    $GENOME_NUM = $clusters{$id}->{num_of_members};

                }

            }

            #getting the next row
            $in = <PANACEA_FF>;

        }

    }

}

#Writting the Javascript for regions (Core and FGR)
sub fgi_javascript() {

    open( OUT, ">", $out_dir . "/$output_id/scripts/fgi.functions.js" );
    my $sz    = $gene_height * 1.0;
    my $sz_h  = $sz;
    my $sz_h2 = $sz_h * 0.5;

    #This is the disk image for the saving both for the fGI page and the core tree page
    my $trans = sprintf( "translate(%f, %f) scale(%f, %f)", 0, 30, 0.045, 0.045 );
    my $action = "onclick=\\\"saveSVGtree(\\\'svg\\\', id)\\\"";

    #disk image with SVG written on it
    my $disk_image_svg = $disk_image_js;
    $disk_image_svg =~ s/TRANS/$trans/g;
    $disk_image_svg =~ s/ACTION/$action/g;
    $disk_image_svg =~ s/TEXT/SVG/g;
    $disk_image_svg =~ s/FS/250/g;

    $trans = sprintf( "translate(%f, %f) scale(%f, %f)", 60, 30, 0.045, 0.045 );
    $action = "onclick=\\\"saveSVGtree(\\\'png\\\', id)\\\"";

    #disk image with PNG written on it
    my $disk_image_png = $disk_image_js;
    $disk_image_png =~ s/TRANS/$trans/g;
    $disk_image_png =~ s/ACTION/$action/g;
    $disk_image_png =~ s/TEXT/PNG/g;
    $disk_image_png =~ s/FS/250/g;

    print OUT "var curr_scroll_top = 0; //Keeps track of where the y-axis scroll bar top is. Used to connect
		var curr_scroll_left = 0; //Same as above but for horizontal
		var curr_gene; //current genome highlighted on tree
		var curr_type; //whether the tree is circular or rectangular
		var rowHighlight = null; //keeps what
		var saveType = \"svg\"; //tores which graphics type is to be saved
		var num_selected = 0; //number of rows selected by the user
		var tree_level; //The depth of the tree- where 1 is the level of the leaves
		var treeJSON = null; //Tree JSON is either null if there is no phylogeny or an array if there is
		var genome2node; //an array where the genome ID is the key and the value is a colon-seperated list of nodes for which that genome is a child
		var node2genome = []; //an array where a Node ID is the key and the value is a colon-seperated list of all the genomes descending from that node
		var metadata_head; //a list of all the metadata variables included
		var metadata; //an array pointing from the metadata values
		var node2meta; //an array that points from the nodes to the associated values
		var nodeIsClicked = 0; //meta-variable keeping track of whether the node is clicked
		var meta_color = []; //A list of the colors associated with each metadata value
		var metadata_is_checked = []; //keeping track of which metadata variables are turned on and off
		var all_tree_nodes = []; //A list of all the nodes in the tree, their level, and the # in the level
		var geneList; //A list of genomes and their colors
		var fgiStart; //A list of 5' bounding genes by fGI
		var fgiEnd; //A list of 3' bounding genes by fGI
		var fgiCnt; //A list of the number genomes by fGI
		var fgiOn; //array showing whether the fGI is shown or not
		var curr_genome = null; //Keeps track of the current gene or fGI turned on and if none are, it is set as null
		var tree_size = $tree_size; //size in pixels of the window panel to show the tree

		if (typeof cntJSON !== 'undefined')
		{
			fgiCnt = JSON.parse(cntJSON);
		}

		if (typeof startJSON !== 'undefined')
		{
			fgiStart = JSON.parse(startJSON);
		}
		if (typeof endJSON !== 'undefined')
		{
			fgiEnd = JSON.parse(endJSON);
		}
		if (typeof geneListJSON !== \'undefined\')
		{
			geneList = JSON.parse(geneListJSON);
		}


		//If there is a phylogeny associated with the genome, load it into treeJSON
		//Also loads the metadata as available
		if (typeof tree_format !== \'undefined\')
		{
			treeJSON = JSON.parse(tree_format);
			for (var i in treeJSON)
			{
				for (var j in treeJSON[i])
				{
					//If the variable is a node, add to array of all the nodes in the tree with the appriopate level (i -> 0) and count in the level (j -> 1)
					if (treeJSON[i][j][6] == \"2\")
					{
						all_tree_nodes[treeJSON[i][j][7]] = [];
						all_tree_nodes[treeJSON[i][j][7]][0] = i;
						all_tree_nodes[treeJSON[i][j][7]][1] = j;
					}
				}
			}

			//Getting the metadata information from the json strings
			genome2node = JSON.parse(gene2node);
			if (typeof meta_data != \'undefined\')
			{
				metadata_head = JSON.parse(meta_head_var);
				metadata = JSON.parse(meta_data);
				node2meta = JSON.parse(meta_node_data);

				//Turning all the metadata checked buttons off at the begining
				for (var i in metadata_head)
				{
					metadata_is_checked[metadata_head[i]] = false;
				}
			}

			//Get the list of genomes that is child to a given node
			for (var i in genome2node)
			{
				var nodes = genome2node[i].split(\";\");
				for (var j = 0; j < nodes.length; j++)
				{
					if (!nodes[j] in node2genome || typeof node2genome[nodes[j]] == \'undefined\')
					{
						node2genome[nodes[j]] = i;
					}
					else
					{
						node2genome[nodes[j]] = node2genome[nodes[j]] + \",\" + i;

					}
				}
			}
		}

		//Draws the page my starting the tree
		function initPage(level, type)
		{
			draw_tree(level, type);

		}
		
		//Changes the size of the Core Region View
		function transformGenes(id, perc)
		{
			var idObj = document.getElementById(\"div\" + id);
			var ratio = parseFloat(idObj.getAttribute(\"transratio\")) * (parseFloat(perc)/100);
			idObj.setAttribute(\"transratio\", ratio);
			var divht = parseFloat(idObj.getAttribute(\"ht\"));
			var divwt = parseFloat(idObj.getAttribute(\"wt\"));
			
			var idObj = document.getElementById(id);
			idObj.setAttribute(\"width\", ratio * divwt);
			idObj.setAttribute(\"height\", ratio * divht);
			idObj.style.transformOrigin = \"0 0\";
			
			var idObj = document.getElementsByClassName(id);
			if (idObj.length > 0)
			{
				for (var i =0; i < idObj.length; i++)
				{
					
					idObj[i].style.transform = \"scale(\" + ratio + \",\"+ratio+\")\";
				}
			}
		}
		
		//Function to move the end-most nodes shown on the tree
		function changeTree(level)
		{
			//Moves level down ie towards the leaves
			if (level == \"-\" && tree_level > 1)
			{
				tree_level = tree_level -1;
				draw_tree(tree_level, curr_type);

				//Both of the below make sure that any previously added information is not lost when the level is changed
				if (curr_genome != null)
				{
					add_fgi_to_tree(curr_genome);
				}

				if (nodeIsClicked != 0)
				{
					restoreNode(nodeIsClicked);
				}

				var p = document.getElementById(\"LevelCount\");
				p.innerHTML = tree_level;
			}

			//Moves levels up ie towards the root
			if (level == \"+\" && tree_level < Object.keys(treeJSON).length)
			{
				tree_level = tree_level +1;
				draw_tree(tree_level, curr_type);

				//Again both
				if (curr_genome != null)
				{
					add_fgi_to_tree(curr_genome);
				}

				if (nodeIsClicked != 0)
				{
					restoreNode(nodeIsClicked);
				}
				var p = document.getElementById(\"LevelCount\");
				p.innerHTML = tree_level;
			}
		}

		//When the metadata type is clicked run through
		//1. Turn click on or off depending
		//2. Change the legend appropiately
		//3. redraw the tree

		function change_metadata_status(evt)
		{
			var target = evt.target;
			metadata_is_checked[target.value] = target.checked;
			change_legend();
			draw_tree(tree_level, curr_type);

		}

		//Change the metadata legend value
		//Use colored circles to show the metadata values
		function change_legend()
		{
			var newMod = document.getElementById(\"LegendSVG\");
			newMod.innerHTML = \"\";
			var tdMeta = document.getElementById(\"SelectMetatype\");

			//Only shows the checked on metadatas
			for (var i in metadata_is_checked)
			{
				var curr_x = 10;

				if (metadata_is_checked[i])
				{
					var metaCheck = document.getElementById(\"select\"+i);

					//Sets the y-variable (height) to the checkbox height
					var cur_y = metaCheck.offsetTop;

					//Cycles through the metadata header info
					for (var j in metadata_head[i])
					{
						//Draws circle for each one...
						newModMod = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\")
						newModMod.setAttribute(\"cx\", curr_x+5);
						newModMod.setAttribute(\"cy\", cur_y - 28);
						newModMod.setAttribute(\"r\", 4);
						newModMod.setAttribute(\"fill\", metadata_head[i][j]);

						newMod.appendChild(newModMod);

						//Moves down the line horizontally with each value
						curr_x = curr_x + 10

						newModMod = document.createElementNS(\"http://www.w3.org/2000/svg\", \"text\")
						newModMod.setAttribute(\"x\", curr_x+2);
						newModMod.setAttribute(\"y\", cur_y -22 );
						newModMod.innerHTML = j;
						newMod.appendChild(newModMod);

						//This is gets the bounding box for the text to see how long the legend word is to appriopately place the next value
						var bbox = newModMod.getBBox();
						curr_x = curr_x + bbox.width;

					}


				}

			}
		}

		//Highlights an fGI island row
		function highlightRow(Name)
		{
			var row = document.getElementById(\"highlight\" + Name);
			var circ = document.getElementById(\"circle\" + Name);

			//If no highlight already on -> just turn on this row
			if (rowHighlight == null)
			{
				row.setAttribute(\"fill\", \"#FEFE22\");
				circ.setAttribute(\"fill\", \"yellow\");
				row.setAttribute(\"fill-opacity\", \"0.5\");
				rowHighlight = Name;
			}
			else
			{

				//Turns off old highlight
				//OldRow is the box across the full row across that is turned yellow
				//OldCirc is the circle that is turned on and off for the highlighting
				var oldRow = document.getElementById(\"highlight\" + rowHighlight);
				var oldCirc = document.getElementById(\"circle\" + rowHighlight);

				//turning the row on
				oldRow.setAttribute(\"fill\", \"none\");
				oldRow.setAttribute(\"fill-opacity\", \"0\");
				oldCirc.setAttribute(\"fill\", \"white\")

				//If the old highlight is the same as the new, don't re-turn it on
				if (Name == rowHighlight)
				{
					rowHighlight = null;
				}
				else
				{
					//Else turn on the new one
					row.setAttribute(\"fill\", \"#FEFE22\");
					row.setAttribute(\"fill-opacity\", \"0.5\");
					circ.setAttribute(\"fill\", \"yellow\");

					rowHighlight = Name;
				}
			}
		}

		//Save all the Image as PNG/SVG if type not equals one or just the window type is one
		function saveFullSVG(type, name)
		{
			var event = new Event(\'build\');
			var nw = document.createElement(\"a\");
			var fileType = document.getElementById(\"saveType\");

			//Get all the DIV parts
			var divLeft = document.getElementById(\"leftDiv\");
			var divRight = document.getElementById(\"rightDiv\");
			var divLeg = document.getElementById(\"legendDiv\");
			var divMain = document.getElementById(\"mainDiv\");

			//Get the size values of the different parts of the DIVs
			var totWidth = parseFloat(divMain.getAttribute(\"tot_width\"));
			var totHeight = parseFloat(divMain.getAttribute(\"tot_height\"));
			var leftWidth = parseFloat(divLeft.style.width);
			var legWidth = parseFloat(divLeg.style.width);
			var legHeight = parseFloat(divLeg.style.height);
			var rightWidth = parseFloat(divRight.style.width);

			//Get the SVG Parts
			var svgLeft = document.getElementById(\"leftSVG\");
			var svgRight = document.getElementById(\"rightSVG\");
			var svgLeg = document.getElementById(\"legendSVG\");
			var svgMain = document.getElementById(\"mainSVG\");
			var svgTop = document.getElementById(\"topRow\");
			var svgBottom = document.getElementById(\"bottomRow\");

			//Get the old SVG sizes from the parts to restore after saving
			var old_ht_left = parseFloat(svgLeft.getAttribute(\"height\"));
			var old_ht_right = parseFloat(svgRight.getAttribute(\"height\"));
			var old_ht_main = parseFloat(svgMain.style.height);
			var old_ht_top = parseFloat(svgTop.style.height);
			var old_ht_bottom = parseFloat(svgBottom.style.height);
			var old_ht_leg = parseFloat(svgLeg.getAttribute(\"height\"));
			var old_wd_left = parseFloat(svgLeft.getAttribute(\"width\"));
			var old_wd_right = parseFloat(svgRight.getAttribute(\"width\"));
			var old_wd_main = parseFloat(svgMain.style.width);
			var old_wd_leg = parseFloat(svgLeg.getAttribute(\"width\"));
			var old_wd_top = parseFloat(svgTop.style.width);
			var old_wd_bottom = parseFloat(svgBottom.style.width);


			//if type is 1, this only saves the window; otherwise it saves all the fgi
			if (type == 1)
			{
				var windWidth = window.innerWidth;
				var windHeight = window.innerHeight;

				if (totWidth > windWidth)
				{
					totWidth = windWidth;
				}
				if (totHeight > windHeight)
				{
					totHeight = windHeight;

				}
			}
			else
			{
				totWidth = old_wd_main + leftWidth + rightWidth;
				totHeight = old_ht_main + legHeight  ;
			}


			var topHeight = parseFloat(svgTop.style.height);
			var bottomHeight = parseFloat(svgBottom.style.height);

			//Setting the Left border area
			svgLeft.setAttribute(\"y\", \"0\");
			svgLeft.setAttribute(\"x\", \"0\");
			svgLeft.setAttribute(\"width\",  leftWidth);
			svgLeft.setAttribute(\"height\",  (totHeight - legHeight));

			//Setting the Main SVG area
			svgMain.setAttribute(\"y\", topHeight);
			svgMain.setAttribute(\"x\",leftWidth);
			svgMain.setAttribute(\"width\",totWidth-rightWidth-leftWidth);
			svgMain.setAttribute(\"height\", (totHeight - legHeight-topHeight -bottomHeight));

			//Setting the top SVG area
			svgTop.setAttribute(\"y\", \"0\");
			svgTop.setAttribute(\"x\",leftWidth);
			svgTop.setAttribute(\"width\",totWidth-rightWidth-leftWidth);
			svgTop.setAttribute(\"height\", (topHeight));

			//Setting the bottom SVG area
			svgBottom.setAttribute(\"y\", totHeight-bottomHeight-legHeight);
			svgBottom.setAttribute(\"x\", leftWidth);
			svgBottom.setAttribute(\"width\",totWidth-rightWidth-leftWidth);
			svgBottom.setAttribute(\"height\", (bottomHeight));

			//Setting the Right SVG area
			svgRight.setAttribute(\"y\",\"0\");
			svgRight.setAttribute(\"x\",totWidth-rightWidth);
			svgRight.setAttribute(\"width\", rightWidth);
			svgRight.setAttribute(\"height\", totHeight - legHeight);

			//Setting the Legend Box
			svgLeg.setAttribute(\"y\", totHeight-legHeight);
			svgLeg.setAttribute(\"x\", leftWidth);
			svgLeg.setAttribute(\"width\",totWidth-rightWidth-leftWidth);
			svgLeg.setAttribute(\"height\",legHeight);

			var legRect = document.getElementById(\"legRect\");
			var oldWidth = parseFloat(legRect.getAttribute(\"width\"));
			if (oldWidth > totWidth-rightWidth-leftWidth)
			{
				legRect.setAttribute(\"width\", totWidth-rightWidth-leftWidth);
			}

			//Making sure that the different annotations in the legend are shown
			var curX = 10;
			var curY = 10;
			var maxY = 0;
			for (var i = 0; i < divLeg.getAttribute(\"num\"); i++)
			{
				var tmpText = document.getElementById(\"TextID\"+i);
				var bbWidth = tmpText.getBBox().width;
				var bbHeight = tmpText.getBBox().height;
				if (bbHeight > maxY)
				{
					maxY = bbHeight;
				}
				if (curX + bbWidth > legRect.getAttribute(\"width\"))
				{

					curX = 10;
					curY = curY + maxY + 1;
					MaxY = 0;
				}
				tmpText.setAttribute(\"x\", curX);
				tmpText.setAttribute(\"y\", curY);
				curX = curX + bbWidth+10;
			}


			var location = window.location.href;
			//Getting the file name for the saved image
			var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);

			//Saving the image as a png
			if (saveType==\"png\")
			{
				var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"svg\");
				newSVG.setAttribute(\"xmlns:xlink\",\"http://www.w3.org/1999/xlink\");
				newSVG.setAttribute(\"height\", totHeight);
				newSVG.setAttribute(\"width\",totWidth);
				newSVG.innerHTML = divLeft.innerHTML + divRight.innerHTML  + divMain.innerHTML  + divLeg.innerHTML;
				var svgStr = new XMLSerializer().serializeToString(newSVG);
				var newSVGBBox = newSVG.getBBox();
				var canvas = document.createElement(\"canvas\");
				canvas.width = totWidth * $DPI;
				canvas.height = totHeight * $DPI;
				var ctx = canvas.getContext(\"2d\");
				ctx.scale($DPI, $DPI);
				var img = new Image();
				img.src = \"data:image/svg+xml;base64,\" + window.btoa(svgStr);
				img.onload = function() {
					ctx.drawImage(img, 0, 0);

					nw.setAttribute(\"download\",name + \"$filenamePNG\");
					nw.setAttribute(\"href\", canvas.toDataURL(\"image/png\"));
					nw.setAttribute(\"target\", \"_blank\");

					alert(\"Saved PNG to \" +  + name + \"$filenamePNG\");
					document.body.appendChild(nw); 
					nw.click(); 
					document.body.removeChild(nw);
				}
			}
			else
			{

				//Saving image as SVG
				var svgStr = \"<svg xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" height=\\\"\"+totHeight+\"\\\" width=\\\"\"+totWidth+\"\\\">\"+divLeft.innerHTML + divRight.innerHTML  + divMain.innerHTML  + divLeg.innerHTML +\"</svg>\";
				var inSvg = new Blob([svgStr], {type: \'image/svg+xml\'});

				nw.setAttribute(\"download\", name + \"$filenameSVG\");
				var svgURL = URL.createObjectURL(inSvg);
				nw.setAttribute(\"href\", svgURL);
				nw.setAttribute(\"target\", \"_blank\");
				document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
				alert(\"Saved SVG to \" + name + \"$filenameSVG\");


			}

			//Resetting the image size variables
			var totWidth = parseFloat(divMain.getAttribute(\"tot_width\"));
			var totHeight = parseFloat(divMain.getAttribute(\"tot_height\"));
			svgLeft.setAttribute(\"height\", old_ht_left);
			svgRight.setAttribute(\"height\", old_ht_right);
			svgMain.style.height=  old_ht_main;
			svgLeg.setAttribute(\"height\", old_ht_leg);
			svgLeft.setAttribute(\"width\", old_wd_left);
			svgRight.setAttribute(\"width\", old_wd_right);
			svgMain.style.width= old_wd_main;
			svgLeg.setAttribute(\"width\", old_wd_leg);

			//Remakes the image with the old values
			makeImage();
		}

		//Changes the view in the main pannel and in the boundary region so that the fGI islands and bounding core genes match up
		function changeWithScroll()
		{
			var svgMain = document.getElementById(\"mainSVG\");
			var num = parseInt(svgMain.getAttribute(\"numrow\"));

			var divLeft = document.getElementById(\"left_border\");
			var divRight = document.getElementById(\"right_border\");
			var divLeg = document.getElementById(\"legendDiv\");
			var divMain = document.getElementById(\"mainDiv\");


			var top = divMain.scrollTop;
			//this is the vertical difference
			var diff = top - curr_scroll_top;
			var left = divMain.scrollLeft;

			//Getting the horizontal difference
			var diff_x = left - curr_scroll_left;

			curr_scroll_top = top;
			curr_scroll_left = left;
			var newHeight = divMain.style.height;
			newHeight = top + newHeight;
			//Moving the core bondaries up and down with diff
			divLeft.setAttribute(\"transform\", \"translate(0,\"+-1*top+\")\");
			divRight.setAttribute(\"transform\", \"translate(0,\"+-1*top+\")\");
			if (diff != 0)
			{
				for (var i = 0; i < num; i++)
				{
					var curRow = document.getElementById(\"mainRow\"+i);
					var rowId = curRow.getAttribute(\"rowid\");
					var curHigh = document.getElementById(\"highlight\"+rowId);
					if (top > parseFloat(curHigh.getAttribute(\"y\"))&& (top + $sz_h)<  (parseFloat(curHigh.getAttribute(\"y\")) + parseFloat(curHigh.getAttribute(\"height\"))))
					{
						curRow.setAttribute(\"y\", top);

					}
				}
			}

			//Moving the Gene and FGI information above and below horizontally with main SVG page
			if (diff_x != 0)
			{
				document.getElementById(\"headerText\").setAttribute(\"x\", parseFloat(document.getElementById(\"headerText\").getAttribute(\"x\"))+diff_x);
				document.getElementById(\"footerText\").setAttribute(\"x\", parseFloat(document.getElementById(\"footerText\").getAttribute(\"x\"))+diff_x);

			}
		}

		//When an fRG is selected, all the genomes associated with those fGR are added to the tree
		function add_fgi_to_tree(fgiID)
		{
			var genomeList = JSON.parse(genomeJSON);

			var new_list = \"\";

			//adding genomes to string by adding name followed by a comma
			for (var i in genomeList[fgiID])
			{
				new_list = new_list + genomeList[fgiID][i] + \",\";
			}
			//Trmming the final comma off the list
			new_list.slice(0,-1);

			//Drawing trees with fGI nodes highlighted
			show_nodes_on_tree(new_list);
			curr_fgi = fgiID;
		}

		//Highlighting the nodes on the phylogeny with descendent genomes containing a selected fGI
		function showFGIonTree(evt)
		{
			var target = evt.target;

			if (! fgiOn)
			{
				var fgiID = target.getAttribute(\"fgiid\");
				target.setAttribute(\"fill\", \"yellow\");
				target.setAttribute(\"fill-opacity\", \"0.5\");
				fgiOn = target;
				add_fgi_to_tree(fgiID);
				curr_genome = fgiID;
			}
			else
			{
				if (fgiOn === target)
				{
					target.setAttribute(\"fill\", \"white\");
					target.setAttribute(\"fill-opacity\", \"0.1\");
					fgiOn = null;
					curr_fgi = null;
					draw_tree(tree_level, curr_type);
					curr_genome = null;
				}
				else
				{

					fgiOn.setAttribute(\"fill\", \"white\");
					fgiOn.setAttribute(\"fill-opacity\", \"0.1\");
					target.setAttribute(\"fill\", \"yellow\");
					target.setAttribute(\"fill-opacity\", \"0.5\");
					fgiOn = target;
					var fgiID = target.getAttribute(\"fgiid\");
					curr_genome = null;
					add_fgi_to_tree(fgiID);
					curr_genome = fgiID;

				}
			}

		}

		//From the FGI page -> creates a page with (a) selected fGIs and (b) phylogeny which ; and it generates a new tab with the page
		function loadTreePage()
		{

			//Only can load tree with fGIs
			if (num_selected==0)
			{
				alert(\"No FGIs selected...<br><br>Please Select at least one<br><br>\");
				return(\"0\");
			}

			//number of fGI rows in the Main image
			var num = document.getElementById(\"mainSVG\").getAttribute(\"numrow\");

			//Array showing which fGI rows are to be shown
			rowOn = new Array();
			//Array showing which gene columns are to be shown
			colOn = new Array();
			//Array showing which fGI rows are associated with which gene columns
			colList = new Array();

			//Counting variable for the number of rows
			var cnt = 0;
			var genomeList = JSON.parse(genomeJSON);

			//string used to stores genome information in JSON
			var newGenomesJSON = \"var genomesJSON='{\";

			//Height of the Gene Bars in the FGI summary figure at the bottom of the page
			var barHeight = 14;
			//String with the SVG
			var svg = \"\";

			var cur_x = 5;
			var max_x = 0;
			for (var i = 0; i < num; i++)
			{

				var rowI = document.getElementById(\"mainRow\" +i);
				var rowListI = rowI.getAttribute(\"rowid\");
				var boxI = document.getElementById(\"genomeBox\" + rowListI);

				if (boxI.getAttribute(\"ison\") ==\"On\")
				{
					//Adding text to the fGI image
					svg = svg + \"<text x =\\\"\" + cur_x+ \"\\\" y=\\\"\"+ ((0.5+cnt)*barHeight)+\"\\\" font-size=\\\"12\\\" dominant-baseline=\\\"middle\\\">\" +fgiStart[rowListI] +\"</text>\";
					if (max_x < fgiStart[rowListI].length)
					{
						max_x = fgiStart[rowListI].length;
					}
					//Setting the rowOn array as one
					rowOn[i] = 1;
					cnt = cnt + 1;

					//Adding the genome to the JSON list
					for (var j in genomeList[rowListI])
					{
						newGenomesJSON = newGenomesJSON + \"\\\"\" + genomeList[rowListI][j] + \"\\\",\";
					}
					newGenomesJSON = newGenomesJSON.slice(0, -1);
					newGenomesJSON = newGenomesJSON + \"],\";
				}
				else
				{
					rowOn[i] = 0;
				}

				//If the fGI row is to be shown, adding the row to the appropiate gene column list
				if (rowOn[i] == 1)
				{
					var genes = rowListI.split(\":\");
					var cnt2 = cnt -1;
					for (var j = 0; j < genes.length; j++)
					{
						//Setting the gene column as on
						colOn[genes[j]] = 1;
						//If to account for uninitiated gene column
						if (typeof colList[genes[j]] === \"undefined\")
						{
							colList[genes[j]] = cnt2 + \":\";
						}
						else
						{
							colList[genes[j]] = colList[genes[j]]+ cnt2 + \":\";
						}
					}
				}
			}

			//Limits number of fGIs added per phylogeny page
			if (cnt>20)
			{
				alert(\"Too many FGIs selected...<br><br>Please select at most six<br><br>\");
				return(\"0\");
			}

			var cur_x = cur_x + max_x * 9;

			//Drawing the gene boxes in the fGI image SVG
			for (var i in geneList)
			{
				var geneID = i.substring(3)
				if (geneID in colOn && colOn[geneID] == 1)
				{
					var colsOld = colList[geneID].slice(0,-1).split(\":\");

					for (var j in colsOld)
					{
						var j1 = parseInt(colsOld[j]);
						svg = svg + \"<rect x=\\\"\" + (cur_x+1) + \"\\\" y =\\\"\" + (j1 * barHeight) + \"\\\" width=\\\"6\\\" height=\\\"\"+barHeight+\"\\\" fill=\\\"\"+geneList[i]+\"\\\"/>\";
					}
					cur_x = cur_x + 8;

				}
			}
			cur_x = cur_x + 9;

			cnt = 0;

			//The maximum size of the 5' end to appropiately place the counts
			max_x = 0;

			//Adding the 3' bounding gene to the fGI summary image
			for (var i = 0; i < num; i++)
			{
				var rowI = document.getElementById(\"mainRow\" +i);
				var rowListI = rowI.getAttribute(\"rowid\");
				var boxI = document.getElementById(\"genomeBox\" + rowListI);

				//Only do so if the fGI row is on
				if (boxI.getAttribute(\"ison\") ==\"On\")
				{
					svg = svg + \"<text x =\\\"\" + cur_x+ \"\\\" y=\\\"\"+ ((0.5+cnt)*barHeight)+\"\\\" font-size=\\\"12\\\">\" +fgiEnd[rowListI] +\"</text>\";
					if (max_x < fgiEnd[rowListI].length)
					{
						max_x = fgiEnd[rowListI].length;
					}
					cnt = cnt + 1;
				}
			}
			cur_x = cur_x + max_x * 9;
			cnt = 0;

			//Adding the number of genomes to the fGI summary image
			for (var i = 0; i < num; i++)
			{
				var rowI = document.getElementById(\"mainRow\" +i);
				var rowListI = rowI.getAttribute(\"rowid\");
				var boxI = document.getElementById(\"genomeBox\" + rowListI);

				if (rowOn[i] ==1)
				{
					svg = svg + \"<text x =\\\"\" + cur_x+ \"\\\" y=\\\"\"+ ((0.5+cnt)*barHeight)+\"\\\" font-size=\\\"12\\\">\" +fgiCnt[rowListI] +\"</text>\";
					cnt = cnt + 1;
				}

			}

			cur_x = cur_x + 25;
			cnt = 0;

			//Adding button to click to (a) show the leaves/nodes on the tree with the selected fGI and (b) to highlight the row
			for (var i = 0; i < num; i++)
			{
				var rowI = document.getElementById(\"mainRow\" +i);
				var rowListI = rowI.getAttribute(\"rowid\");
				var boxI = document.getElementById(\"genomeBox\" + rowListI);
				if (rowOn[i] ==1)
				{
					svg = svg + \"<rect x =\\\"0\\\" y=\\\"\"+ ((cnt)*barHeight)+\"\\\" height=\\\"\"+barHeight+\"\\\" width=\\\"\" + cur_x+ \"\\\" fill=\\\"white\\\" fill-opacity=\\\"0.05\\\" id=\\\"rect\" + cnt + \"\\\" fgiID=\\\"\" + rowListI +\"\\\" onclick=\\\"showFGIonTree(evt)\\\"/>\";
					cnt = cnt + 1;
				}
			}

			newGenomesJSON = newGenomesJSON.slice(0,-1)  + \"}';\";

			//Writing the tree html, then opening it in a new page
			//Not saved in any file


			var html = \"<html><script type=\\\"text/javascript\\\">\";

			//Adding JSON string with genomes
			html = html + \"var genomeJSON= \'\" + genomeJSON + \"\';\";
			html = html + \"</script>\";

			//Adding JSON string with the phylogeny written in
			html = html + \"<script type='text/javascript' src='../json/tree.json'></script><script type='text/javascript' src='../scripts/fgi.functions.js'></script>\";
			html = html + \"<body style='font-family:Courier New' xmlns=\\\"http://www.w3.org/1999/xhtml\\\"  xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" oninit=\\\"draw_tree(4, \'circular\')\\\">\";
			html = html + \"<div id='select'>\";

			html = html + \"<table ><tr height=\\\"\" + (tree_size * 0.2) + \"\\\"><td style=\\\"vertical-align:top;\\\" rowspan=\\\"3\\\">\"

			html = html + \"<div><svg version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" id=\\\"tree_svg\\\" height=\\\"\"+ tree_size+\"\\\" width=\\\"\"+ tree_size+\"\\\"></svg>\";

			html = html + \"</div></td><td style=\\\"vertical-align:top;\\\"><div><h3>Tree Functions</h3>Outer Ring Tree Level:<button onclick=\\\"changeTree(\'-\')\\\">-</button><b id = \\\"LevelCount\\\">4</b><button onclick=\\\"changeTree(\'+\')\\\">+</button>\";

			html = html + \"<form>Phylogeny Style:<input type=\\\"radio\\\" name=\\\"type\\\" value=\\\"Circular\\\" checked=\\\"true\\\" onclick=\\\"changeType(event)\\\"/>Circular<input type=\\\"radio\\\" name=\\\"type\\\" value=\\\"Linear\\\" onclick=\\\"changeType(event)\\\"/>Linear</form>\";

			html = html + \"<svg version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" id=\\\"save_svg\\\" height=\\\"\"+ (tree_size *0.1)+\"\\\" width=\\\"\"+ (tree_size*0.4)+\"\\\"><text x=\\\"2\\\" y=\\\"14\\\">Node Colors:</text><circle cx=\\\"155\\\" r=\\\"4\\\" cy=\\\"8\\\" fill=\\\"black\\\"/><text x=\\\"162\\\" y=\\\"14\\\">Present</text><circle cx=\\\"245\\\" r=\\\"4\\\" cy=\\\"8\\\" fill=\\\"LightGray\\\"/><text x=\\\"252\\\" y=\\\"14\\\">Missing</text>$disk_image_svg $disk_image_png</svg></div></td></tr><tr height=\\\"\"+ (tree_size *0.4)+\"\\\" >\";
			html = html + \"<td style=\\\"vertical-align:top;\\\"><div id=\\\"gene_list\\\" style=\\\"height:\"+ (tree_size*0.4)+\"; overflow-y: auto;\\\"></div></td></tr><tr height=\\\"\"+ (tree_size*0.15)+\"\\\"><td style=\\\"vertical-align:bottom;\\\">\";
			html = html + \"<div id=\\\"percentage_image\\\"></div></td></tr><tr height=\\\"\"+ (tree_size*0.25)+\"\\\">\";
			html = html + \"<td><table><tr><td id=\\\"SelectMetatype\\\" width=\\\"\"+ (tree_size*0.4)+\"\\\" style=\\\"vertical-align:top;\\\"><h3>Metadata Type</h3>\";
			for (var a in metadata_head)
			{
				html = html + \"<input type=\\\"checkbox\\\" id=\\\"select\"+a+\"\\\" value=\\\"\" + a+ \"\\\" onchange=\\\"change_metadata_status(event)\\\" />\"+a+\"<p>\";
			}
			html = html +\"</td><td id = \\\"Legend\\\" style=\\\"vertical-align:top;\\\"><h3>Legend</h3><svg id=\\\"LegendSVG\\\" version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" height=\\\"\"+ (tree_size*0.2)+\"\\\" width=\\\"\"+ (tree_size*0.6)+\"\\\"></svg></td></tr></table><td style=\\\"vertical-align:top;\\\"><h3>Percentages</h3><div id =\\\"percDiv\\\"></div></td></tr></table>\";

			html = html + \"<svg id = \\\"selectFGR\\\" version=\\\"1.2\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" y=\\\"0\\\" x=\\\"0\\\" width=\\\"\" + cur_x+ \"\\\" height=\\\"\" + ((cnt) * barHeight)+ \"\\\">\";

			html = html + svg + \"<\\svg>\";

			html = html + \"</div></body>\";
			html = html + \"<script type=\\\"text/javascript\\\">\";
			html = html + \"draw_tree(4, \'circular\'); \";
			html = html + \"</script></html>\";
			var newWind = window.open();
			newWind.document.write(html);
			newWind.document.close();
		}

		//Given a selected set of fGIs, this function removes all other rows of fGIs not selected and all columns of genes not contained in the selected fGIs
		function trimGenomesImage()
		{
			//This changes the title of the trim button to \'reset\' if trimming done
			//Otherwise, undoes trimming (ie resets) and redraws
			var trimBut = document.getElementById(\"trimButton\");

			//If trimming is selected....
			if (trimBut.getAttribute(\"trim\") == \"0\")
			{

				if (num_selected==0)
				{
					alert(\"No FGIs or Genomes selected...<br><br>Please Select at least one\<br><br>\");
					return(\"0\");
				}

				trimBut.setAttribute(\"trim\",\"1\");

				//Renaming the trim button to allow resetting of viewing all of the fGIs
				var trimText = document.getElementById(\"trimButton\");
				trimText.innerHTML = \"Reset\";


				//Array showing which fGI rows are to be shown
				rowOn = new Array();

				//Array showing which gene columns are to be shown
				colOn = new Array();

				//Array matching gene column to a given x-axis value
				colX = new Array();

				//num is the number fGR rows
				var num = document.getElementById(\"mainSVG\").getAttribute(\"numrow\");

				var genomeList = JSON.parse(genomeJSON);

				//This is the location of the top Y-value
				var topAdd = parseFloat(document.getElementById(\"mainSVG\").style.top);

				//Current x and y-axis locations
				var curY = 0;
				var curX = 0;

				for (var i = 0; i < num; i++)
				{
					//Going through the rows, getting the add the fGR row information
					var rowI = document.getElementById(\"mainRow\" +i);
					var rowListI = rowI.getAttribute(\"rowid\");
					var rowLi = document.getElementById(\"rowLeft\" + i)

					//Also look to any shown fGRS individual genomes
					var boxI = document.getElementById(\"genomeBox\" + rowListI);
					var showI = document.getElementById(\"All\"+rowListI);
					var highI = document.getElementById(\"highlight\" + rowListI);

					//If the fGR box is checked, include it
					if (boxI.getAttribute(\"ison\") ==\"On\")
					{
						rowOn[i] = 1;
					}
					else
					{
						rowOn[i] = 0;

						//including selected genomes as well
						if (showI.getAttribute(\"ison\") == \"On\")
						{
							var genomes = genomeList[rowListI];
							for (var j = 0; j < genomes.length; j++)
							{
								var boxGenomeI = document.getElementById(\"genomeBox\" + genomes[j]);
								if (boxGenomeI.getAttribute(\"ison\") == \"On\")
								{
									rowOn[i] = 1;
								}
							}
						}
					}

					//If the fGR row is turned on, start counting the y-axis location and labeling the genes in the fGR as on
					//Resetting the y-axis to the appropiate height for the row and buttons
					if (rowOn[i] == 1)
					{
						var genes = rowListI.split(\":\");
						for (var j = 0; j < genes.length; j++)
						{
							colOn[genes[j]] = 1;
						}
						highI.setAttribute(\"y\", curY);
						showI.setAttribute(\"y\", curY + topAdd);
						rowI.setAttribute(\"y\", curY);
						rowLi.setAttribute(\"y\", curY + topAdd);
						curY = curY + parseFloat(showI.getAttribute(\"height\"));

						//Doing the same for the individual genomes shown
						var genomes = genomeList[rowListI];
						if (showI.getAttribute(\"ison\") == \"On\")
						{
							for (var j = 0; j < genomes.length; j++)
							{

								var svgGenomeI = document.getElementById(\"genomeSVG\" + genomes[j]);
								var boxGenomeI = document.getElementById(\"genomeBox\" + genomes[j]);

								if (boxGenomeI.getAttribute(\"ison\") == \"On\")
								{
									svgGenomeI.setAttribute(\"y\", curY + topAdd)
									curY = curY + parseFloat(svgGenomeI.getAttribute(\"height\"));
								}
								else
								{
									//Removing all the genomes not turned on, by setting height to zero
									svgGenomeI.setAttribute(\"height\", 0)
								}

							}
						}
					}
					else
					{
						//Removing all the FGRs not turned on, by setting height to zero
						showI.setAttribute(\"height\", 0);
						rowI.setAttribute(\"height\", 0);
						rowLi.setAttribute(\"height\", 0);
						highI.setAttribute(\"height\", 0);
					}
				}

				//SVG elements that need to have their heights and location reset
				var mainSVG = document.getElementById(\"mainSVG\");
				var topRow = document.getElementById(\"topRow\");
				var botRow = document.getElementById(\"bottomRow\");


				mainSVG.style.height =  curY;
				botRow.style.top = curY+topAdd;
				var group = botRow.getElementsByTagName(\"rect\");
				curX = 0;

				//Goes through the genes to keep the gene names and region that are kept in the top and bottom row
				//Also keeps the arrow images in the main image if the column and fGI is turned on
				if (group != null)
				{
					for (var i = 0; i < group.length; i++)
					{
						id = group[i];
						if (id.hasAttribute(\"id\"))
						{
							var geneNum = id.getAttribute(\"id\").split(\"_\");
							var topRect = document.getElementById(\"topCL_\" + geneNum[1]);
							var mainRect = document.getElementById(\"mainCL_\" + geneNum[1]);
							var topText = document.getElementById(\"topTextCL_\" + geneNum[1]);
							var botText = document.getElementById(\"botTextCL_\" + geneNum[1]);
							var lineID = document.getElementById(\"lineCL_\" + geneNum[1]);
							if (colOn[geneNum[1]] == 1)
							{
								id.setAttribute(\"oldX\", id.getAttribute(\"x\"));
								topRect.setAttribute(\"oldX\", topRect.getAttribute(\"x\"));
								mainRect.setAttribute(\"oldX\", mainRect.getAttribute(\"x\"));
								id.setAttribute(\"x\", curX);
								topRect.setAttribute(\"x\", curX);
								mainRect.setAttribute(\"x\", curX);
								var newCurX =  curX + 12.5;

								//Storing the old values
								topText.setAttribute(\"oldX\",topText.getAttribute(\"x\"));

								topText.setAttribute(\"oldTransform\", topText.getAttribute(\"transform\"));
								botText.setAttribute(\"oldX\",botText.getAttribute(\"x\"));
								botText.setAttribute(\"oldTransform\", botText.getAttribute(\"transform\"));

								//setting the new x-axis values for the top and bottom IDs and the line
								topText.setAttribute(\"x\",newCurX);

								topText.setAttribute(\"transform\", \"rotate(270,\" +newCurX + \",\"+ topText.getAttribute(\"y\")+\")\");
								botText.setAttribute(\"x\",newCurX);
								botText.setAttribute(\"transform\", \"rotate(-90,\" + newCurX+ \",\"+ botText.getAttribute(\"y\")+\")\");

								lineID.setAttribute(\"x1\", newCurX);
								lineID.setAttribute(\"x2\", newCurX);

								colX[geneNum[1]] = curX;
								curX = curX + parseFloat(id.getAttribute(\"width\"));
							}
							else
							{
								//if turned off, set all widths to zero
								//First through keeping all the old values
								id.setAttribute(\"oldWidth\", id.getAttribute(\"width\"));
								topRect.setAttribute(\"oldWidth\", topRect.getAttribute(\"width\"));
								mainRect.setAttribute(\"oldWidth\", mainRect.getAttribute(\"width\"));
								id.setAttribute(\"oldX\", id.getAttribute(\"x\"));
								topRect.setAttribute(\"oldX\", topRect.getAttribute(\"x\"));
								mainRect.setAttribute(\"oldX\", mainRect.getAttribute(\"x\"));

								//know setting the widths to zero
								id.setAttribute(\"width\", 0);
								topRect.setAttribute(\"width\", 0);
								mainRect.setAttribute(\"width\", 0);
								id.setAttribute(\"x\", 0);
								topRect.setAttribute(\"x\", 0);
								mainRect.setAttribute(\"x\", 0);
								topText.setAttribute(\"oldText\", topText.innerHTML);
								botText.setAttribute(\"oldText\", botText.innerHTML);
								topText.innerHTML = \"\";
								botText.innerHTML = \"\";
								lineID.setAttribute(\"stroke\", \"none\");
							}
						}
					}
				}
				//Array of genes with information
				var geneList = JSON.parse(fastaJSON);

				//Resetting the x-axis location of each gene arrow in the main image
				for (var i = 0; i < num; i++)
				{
					var rowI = document.getElementById(\"mainRow\" +i);
					var rowListI = rowI.getAttribute(\"rowid\");
					var rowLi = document.getElementById(\"rowLeft\" + i)
					var boxI = document.getElementById(\"genomeBox\" + rowListI);
					var showI = document.getElementById(\"All\"+rowListI);

					if (rowOn[i] == 1)
					{
						var geneList = rowListI.split(\":\");
						for (var j = 0; j < geneList.length; j++)
						{
							var tmpPath = document.getElementById(i + \"_CL_\" + geneList[j]);
							if (tmpPath != null)
							{
								tmpPath.setAttribute(\"x\", colX[geneList[j]]);
							}
						}
					}
				}

				//Resetting document sizes and locations before remaking
				document.getElementById(\"headerText\").setAttribute(\"x\", curX/2);
				document.getElementById(\"footerText\").setAttribute(\"x\", curX/2);
				mainSVG.style.height = curY;
				mainSVG.style.width = curX;
				topRow.style.width = curX;
				botRow.style.width = curX;

				var mainDiv = document.getElementById(\"mainDiv\");
				var totHeight = curY + parseFloat(topRow.style.height) + parseFloat(botRow.style.height);

				mainDiv.style.height = totHeight;
				mainDiv.style.width = curX;
				mainDiv.setAttribute(\"tot_width\", parseFloat(mainDiv.style.width)+parseFloat(document.getElementById(\"leftDiv\").style.width)+parseFloat(document.getElementById(\"rightDiv\").style.width));
				mainDiv.setAttribute(\"tot_height\", totHeight+parseFloat(document.getElementById(\"legendDiv\").style.height));

				document.getElementById(\"legendDiv\").style.top = totHeight;
				document.getElementById(\"leftDiv\").style.height = totHeight;
				document.getElementById(\"rightDiv\").style.height = totHeight;
				document.getElementById(\"saveDiv\").style.top = totHeight;
				document.getElementById(\"sortDiv\").style.top = totHeight;

				//Remaking the image
				makeImage();
			}
			else
			{
				//Resetting the image
				location.reload();
			}

		}

		//Picks the genomes once the
		function isFastaOn(name, type)
		{
			var box = document.getElementById(\"genomeBox\" + name);
			var text = document.getElementById(\"genomeTextBox\" + name);

			if (type == \"Full\")
			{
				//For fGR regions
				var genomeList = JSON.parse(genomeJSON);

				var genomes = genomeList[name];
				if (box.getAttribute(\"ison\") == \"Off\")
				{
					text.setAttribute(\"fill\", \"black\"); //adding the x in the box to show that it is turned on
					box.setAttribute(\"ison\",\"On\");
					for (var i = 0; i < genomes.length; i++)
					{
						if (isGenomeOn[genomes[i]] == 0)
						{
							num_selected = num_selected + 1;
						}
						isGenomeOn[genomes[i]] = 1;

						//Turning all the genomes in the fGR on as well
						var box2 = document.getElementById(\"genomeBox\" + genomes[i]);
						var textBox2 = document.getElementById(\"genomeTextBox\" + genomes[i]);

						if (box2 != null)
						{
							box2.setAttribute(\"ison\", \"On\");
							textBox2.setAttribute(\"fill\", \"black\");
						}
					}
				}
				else
				{
					text.setAttribute(\"fill\", \"white\"); //removing the x in the box to show that it is turned off
					box.setAttribute(\"ison\",\"Off\");
					for (var i = 0; i < genomes.length; i++)
					{
						if (isGenomeOn[genomes[i]] == 1)
						{
							num_selected = num_selected - 1;
						}
						isGenomeOn[genomes[i]] = 0;

						//Turning all the genomes in the fGR off as well

						var box2 = document.getElementById(\"genomeBox\" + genomes[i]);
						var textBox2 = document.getElementById(\"genomeTextBox\" + genomes[i]);

						if (box2 != null)
						{
							box2.setAttribute(\"ison\", \"Off\");
							textBox2.setAttribute(\"fill\", \"white\");
						}
					}
				}
			}
			else
			{
				//Single genomes: same on and off
				if (isGenomeOn[name] == 0)
				{
					isGenomeOn[name] = 1;
					num_selected = num_selected + 1;

					var box2 = document.getElementById(\"genomeBox\" +name);
					var textBox2 = document.getElementById(\"genomeTextBox\" + name);
					if (box2 != null)
					{
						box2.setAttribute(\"ison\", \"On\");
						textBox2.setAttribute(\"fill\", \"black\");
					}
				}
				else
				{
					isGenomeOn[name] = 0;
					num_selected = num_selected - 1;

					var box2 = document.getElementById(\"genomeBox\" +name);
					var textBox2 = document.getElementById(\"genomeTextBox\" + name);
					if (box2 != null)
					{
						box2.setAttribute(\"ison\", \"Off\");
						textBox2.setAttribute(\"fill\", \"white\");
					}
				}
			}

		}

		//Shows the genomes associated with a given fGI
		function showGenomes(Name)
		{
			//SVG names
			var nameSVG = document.getElementById(\"All\"+Name);
			var svgMain = document.getElementById(\"mainSVG\");
			var isOn = nameSVG.getAttribute(\"ison\");

			var showText = document.getElementById(\"showText\"+Name);

			var svgRight = document.getElementById(\"rightSVG\");
			var svgMain = document.getElementById(\"mainSVG\");
			var svgLeft = document.getElementById(\"leftSVG\");
			var bottowRow = document.getElementById(\"bottomRow\");
			var topRow = document.getElementById(\"topRow\");

			var svgRight = document.getElementById(\"rightSVG\");
			var svgMain = document.getElementById(\"mainSVG\");
			var svgLeft = document.getElementById(\"leftSVG\");
			var bottowRow = document.getElementById(\"bottomRow\");

			//get all the rects (highlight boxes) and lines
			var rects = svgMain.getElementsByTagName(\"rect\");
			var lines = svgMain.getElementsByTagName(\"line\");

			var Names = JSON.parse(namesJSON);
			var num = parseInt(nameSVG.getAttribute(\"numrow\"));
			var rowNum = parseInt(svgMain.getAttribute(\"numrow\"));

			//If the extended genome view is turned off...
			if (isOn == \"Off\")
			{

				var genomeList = JSON.parse(genomeJSON);
				var genomes = genomeList[Name];
				var curY = parseFloat(nameSVG.getAttribute(\"height\"))+$sz_h2;
				var curY2 = $sz_h2;

				//For each of the genomes that has the fGI path (1) make an SVG with genome name and selection box
				for (var i = 0; i < genomes.length; i++)
				{
					var genomeSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"svg\");
					genomeSVG.setAttribute(\"height\",\"$sz_h2\");
					genomeSVG.setAttribute(\"y\",curY-$sz_h2);
					genomeSVG.setAttribute(\"x\",0);
					genomeSVG.setAttribute(\"height\",$sz_h2);

					genomeSVG.setAttribute(\"id\",\"genomeSVG\" + genomes[i]);

					var genomeText = document.createElementNS(\"http://www.w3.org/2000/svg\", \"text\");
					genomeText.innerHTML = genomes[i];
					genomeText.setAttribute(\"height\",\"$sz_h2\");

					genomeText.setAttribute(\"y\",$sz_h2);
					genomeText.setAttribute(\"x\",2 * $sz_h);
					genomeText.setAttribute(\"font-color\",\"black\");
					genomeText.setAttribute(\"font-size\",$sz_h2);
					genomeSVG.setAttribute(\"type\",\"genome\");
					genomeSVG.appendChild(genomeText);

					var genomeBox = document.createElementNS(\"http://www.w3.org/2000/svg\", \"svg\");
					genomeBox.setAttribute(\"onclick\", \"isFastaOn(\'\"+genomes[i]+\"\', \'single\')\");
					genomeBox.setAttribute(\"visibility\", \"visible\");
					genomeBox.setAttribute(\"id\", \"genomeBox\"+genomes[i]);

					genomeBox.setAttribute(\"y\", 0);
					genomeBox.setAttribute(\"x\", (1) * $sz_h);
					genomeBox.setAttribute(\"type\",\"genome\");


					var rectBox = document.createElementNS(\"http://www.w3.org/2000/svg\", \"rect\");
					var textBox = document.createElementNS(\"http://www.w3.org/2000/svg\", \"text\");

					rectBox.setAttribute(\"fill\", \"white\");
					rectBox.setAttribute(\"stroke\", \"black\");
					rectBox.setAttribute(\"y\", $sz_h2 * 0.1);
					rectBox.setAttribute(\"x\", $sz_h2 * 0.1);
					rectBox.setAttribute(\"height\", $sz_h2 * 0.8);
					rectBox.setAttribute(\"width\", $sz_h2 * 0.8);
					rectBox.setAttribute(\"id\", \"genomeTextBox\"+genomes[i]);
					

					/* Below used to be for the x to be checked
					textBox.setAttribute(\"font-size\", $sz_h2 * .75);
					textBox.setAttribute(\"y\", $sz_h2 * 0.5);
					textBox.setAttribute(\"x\", $sz_h2 * 0.5);
					textBox.setAttribute(\"text-anchor\", \"middle\");
					textBox.setAttribute(\"alignment-baseline\", \"central\");
					*/
					//If the genome is \"on\" ie checked- add the x to the text box otherwise add blank space
					if (isGenomeOn[genomes[i]]==1)
					{
						textBox.getAttribute(\"fill\", \"black\");
						genomeBox.setAttribute(\"ison\", \"true\");
					}
					else
					{
						textBox.getAttribute(\"fill\", \"white\");
						genomeBox.setAttribute(\"ison\", \"false\");

					}

					//Add check box to inner SVG and that to the total SVG
					genomeBox.appendChild(rectBox);
					genomeBox.appendChild(textBox);
					genomeSVG.appendChild(genomeBox);
					nameSVG.appendChild(genomeSVG);

					curY = curY + $sz_h2;
					curY2 = curY2 + $sz_h2;
				}
				//Change the text of the choose genomes from + to -
				showText.setAttribute(\"stroke\", \"white\");
				//Also change the highlight box height and add broder to make more visible
				var curHigh = document.getElementById(\"highlight\"+Name);
				curHigh.setAttribute(\"height\", curY2+$sz_h2);
				curHigh.setAttribute(\"stroke\", \"#888888\");


				nameSVG.setAttribute(\"height\", curY);
				nameSVG.setAttribute(\"ison\", \"On\");
				nameSVG.setAttribute(\"addedY\", curY2);
				curY = curY + parseFloat(nameSVG.getAttribute(\"y\"));

				//Going through all the lower rows to add to the appropiate heights
				for (var i = num+1; i < rowNum; i++)
				{
					var nextSVG = document.getElementById(\"All\"+Names[i]);
					nextSVG.setAttribute(\"y\", curY);
					var curRow = document.getElementById(\"mainRow\"+i);
					var curLeft = document.getElementById(\"rowLeft\"+i);
					var curHigh = document.getElementById(\"highlight\"+Names[i]);
					curRow.setAttribute(\"y\", parseFloat(curRow.getAttribute(\"y\"))+curY2);
					curRow.setAttribute(\"old_y\", parseFloat(curRow.getAttribute(\"y\"))+curY2);

					curLeft.setAttribute(\"y\", curY);
					curHigh.setAttribute(\"y\", parseFloat(curHigh.getAttribute(\"y\"))+curY2);
					curY = curY +parseFloat(nextSVG.getAttribute(\"height\"));
				}

				//Doing the same for the highlight rectangle
				for (var i = 0; i < rects.length; i++)
				{
					if(rects[i].getAttribute(\"type\")==\"background\")
					{
						rects[i].setAttribute(\"height\",curY);
					}
				}

				//Extend the vertical lines to the new lengths
				for (var i = 0; i < lines.length; i++)
				{
					lines[i].setAttribute(\"height\",curY);
				}
				svgRight.setAttribute(\"height\", curY);
				svgMain.style.height = curY;
				svgLeft.setAttribute(\"height\", curY);
				bottomRow.style.top = curY;

				//Redraws image
				makeImage();
			}

			//If fGR is already showing all the genomes, close this by removing it all
			if (isOn == \"On\")
			{
				//Finds all the SVGs associated with the and removes them
				var texts = nameSVG.getElementsByTagName(\"svg\");
				for (var i = 0; i < texts.length; i++)
				{
					var Parent = texts[i].parentNode;
					if (texts[i].getAttribute(\"type\") == \"genome\")
					{
						Parent.removeChild(texts[i]);
						i = i - 1;
					}
				}

				//Gets the added Y to remove from all the lower fgis
				var curY = nameSVG.getAttribute(\"addedY\");

				nameSVG.removeAttribute(\"addedY\");
				var curHigh = document.getElementById(\"highlight\"+Name);
				curHigh.setAttribute(\"height\", $sz_h);
				curHigh.setAttribute(\"stroke\", \"none\");

				numRow = document.getElementById(\"mainRow\" + curHigh.getAttribute(\"num\"));
				numRow.setAttribute(\"y\", parseFloat(curHigh.getAttribute(\"y\")));


				for (var i = num+1; i < rowNum; i++)
				{
					var nextSVG = document.getElementById(\"All\"+Names[i]);
					var curRow = document.getElementById(\"mainRow\"+i);
					var curLeft = document.getElementById(\"rowLeft\"+i);
					var curHigh = document.getElementById(\"highlight\"+Names[i]);
					nextSVG.setAttribute(\"y\",parseFloat(nextSVG.getAttribute(\"y\"))-curY);

					curRow.setAttribute(\"y\", parseFloat(curRow.getAttribute(\"y\"))-curY);
					curLeft.setAttribute(\"y\", curLeft.getAttribute(\"y\") - curY);
					curHigh.setAttribute(\"y\", parseFloat(curHigh.getAttribute(\"y\"))-curY);
				}

				//Resetting the overall sizes
				svgRight.setAttribute(\"height\", parseFloat(svgRight.getAttribute(\"height\")) - curY);
				svgMain.style.height = parseFloat(svgMain.style.height) - curY;
				svgLeft.setAttribute(\"height\", parseFloat(svgLeft.getAttribute(\"height\")) - curY);
				bottomRow.style.top = parseFloat(bottomRow.style.top) - curY;

				for (var i = 0; i < rects.length; i++)
				{
					if(rects[i].getAttribute(\"type\")==\"background\")
					{
						rects[i].setAttribute(\"height\",parseFloat(rects[i].getAttribute(\"height\")) - curY);
					}
				}
				for (var i = 0; i < lines.length; i++)
				{
					lines[i].setAttribute(\"height\",parseFloat(lines[i].getAttribute(\"height\")) - curY);
				}
				showText.setAttribute(\"stroke\", \"black\");
				nameSVG.setAttribute(\"height\", parseFloat(nameSVG.getAttribute(\"height\"))- curY);
				nameSVG.setAttribute(\"ison\", \"Off\");

				//Redrawing image
				makeImage();

			}

		}

		//Redraws the tree
		function changeType(event)
		{
			var target = event.target;
			draw_tree(tree_level, target.value);
		}

		//Turns the phylogeny node on or off
		function keepNodeOn(gene_id)
		{
			if (nodeIsClicked == 0)
			{
				nodeIsClicked = gene_id;
			}
			else
			{
				nodeIsClicked = 0;
			}
		}


		function score_tree(gene_id)
		{
			if (curr_gene != null)
			{
				var old_curr_gene = curr_gene;
				curr_gene = null;
				var old_tar = document.getElementById(old_curr_gene + \"detailHT\");
				old_tar.setAttribute(\"stroke\", \"black\");
				old_tar.setAttribute(\"stroke-width\", \"1\")
				var old_tar = document.getElementById(old_curr_gene + \"detail\");
				old_tar.setAttribute(\"stroke\", \"black\");
				old_tar.setAttribute(\"stroke-width\", \"1\")

				draw_tree(tree_level, curr_type);
				if (old_curr_gene == gene_id)
				{
					return; //exits the function
				}
			}
			curr_gene = gene_id;
			var old_tar = document.getElementById(gene_id + \"detailHT\");
			old_tar.setAttribute(\"stroke\", \"black\");
			old_tar.setAttribute(\"stroke-width\", \"3\")
			var old_tar = document.getElementById(gene_id + \"detail\");
			old_tar.setAttribute(\"stroke\", \"black\");
			old_tar.setAttribute(\"stroke-width\", \"3\")

			add_to_tree(gene_id);
		}


		function add_to_tree(gene_id)
		{
			var fasta = JSON.parse(fastaJSON);
			var treeSVG = document.getElementById(\"tree_svg\");

			curr_gene = gene_id;

			if (gene_id in fasta)
			{
				var genome_list = fasta[gene_id][\"Genomes\"];
				var new_list = genome_list.replace(/<br>/g, \"\");
				show_nodes_on_tree(new_list);

			}
		}

		function show_nodes_on_tree(new_list)
		{
			var treeSVG = document.getElementById(\"tree_svg\");
			curr_fgi = null;
			var tmp_curr_gene = curr_gene;
			curr_gene = null;
			draw_tree(tree_level, curr_type);

				var genomes = new_list.split(\",\");
				var node_counts = [];
				for (i in genomes)
				{
					if (genomes[i] in genome2node)
					{
						var node_list = genome2node[genomes[i]].split(\";\");
						for (j in node_list)
						{
							if (node_list[j] in node_counts)
							{
								node_counts[node_list[j]] = node_counts[node_list[j]] + 1;
							}
							else
							{
								node_counts[node_list[j]] =1;
							}
						}
					}
				}
				for (j in node2genome)
				{
					if (j in node_counts)
					{
						var line = document.getElementsByClassName(j + \"0\");
						if (line != null)
						{
							for (var k =0 ;k < line.length; k++)
							{
								line[k].setAttributeNS(null,\"stroke-width\", 4);
							}
						}
						var line = document.getElementsByClassName(j + \"1\");
						if (line != null)
						{
							for (var k =0 ;k < line.length; k++)
							{

								line[k].setAttributeNS(null,\"stroke-width\", 4);
							}
						}
						var circ = document.getElementsByClassName(j + \"2\");
						if (circ != null)
						{
							for (var k =0 ;k < circ.length; k++)
							{

							var tot = parseInt(circ[k].getAttributeNS(null, \"cnt\"));
							if (node_counts[j] === tot)
							{
								if (j != nodeIsClicked)
								{
									circ[k].setAttributeNS(null, \"fill\", \"black\");

								}
								var old_r = parseFloat(circ[k].getAttribute(\"r\"));
								circ[k].setAttributeNS(null, \"r\", old_r * 2);
								circ[k].setAttributeNS(null, \"old_r\", old_r);
								circ[k].setAttributeNS(null, \"totCnt\", tot);
								circ[k].setAttributeNS(null, \"genomeCnt\", node_counts[j]);
							}
							else
							{
								circ[k].setAttributeNS(null, \"stroke\", circ[k].getAttribute(\"fill\"));

								circ[k].setAttributeNS(null, \"fill\", \"LightGray\");

							var cx = parseFloat(circ[k].getAttributeNS(null, \"cx\"));
								var cy = parseFloat(circ[k].getAttributeNS(null, \"cy\"));
								var r = parseFloat(circ[k].getAttributeNS(null, \"r\"))*2;
								circ[k].setAttributeNS(null, \"r\", r);
								circ[k].setAttributeNS(null, \"totCnt\", tot);
								circ[k].setAttributeNS(null, \"genomeCnt\", node_counts[j]);
								var largeArc = 0;
								if ((node_counts[j] / tot) > 0.5)
								{
									largeArc = 1;
								}
								var x1 = Math.cos(0) * r + cx;
								var y1 = Math.sin(0) * r + cy;
								var x2 = Math.cos(node_counts[j] / tot * 2 * Math.PI) * r + cx;
								var y2 = Math.sin(node_counts[j] / tot * 2 * Math.PI) * r + cy;

								var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"path\");
								newSVG.setAttributeNS(null, \"id\", j + \"3\");
								newSVG.setAttributeNS(null, \"d\", \"M\" + cx+\" \" +cy+\" L \" + x1+\" \" +y1+\" A \" + r + \" \" + r + \" 0 \"+ largeArc+ \" 1 \" +x2+\" \" + y2 + \" L \" + cx + \" \" + cy);
								newSVG.setAttributeNS(null, \"stroke\", \"red\");
								newSVG.setAttributeNS(null, \"stroke-width\", \"0\");
								newSVG.setAttributeNS(null, \"fill\", \"black\");
								newSVG.setAttributeNS(null, \"onmouseover\", circ[k].getAttributeNS(null, \"onmouseover\"));
								newSVG.setAttributeNS(null, \"onmouseout\",circ[k].getAttributeNS(null, \"onmouseout\"));
								newSVG.setAttributeNS(null, \"onclick\", circ[k].getAttributeNS(null, \"onclick\"));

								if (j == nodeIsClicked)
								{
									newSVG.setAttributeNS(null, \"fill\", \"green\");
								}
								treeSVG.appendChild(newSVG);

							}
							}
						}

					}
					else
					{
						var circ = document.getElementsByClassName(j + \"2\");
						for (var k =0 ;k < circ.length; k++)
						{

							var old_r = parseFloat(circ[k].getAttribute(\"r\"));
							circ[k].setAttributeNS(null, \"r\", old_r * 2);
							circ[k].setAttributeNS(null, \"stroke\", circ[k].getAttribute(\"fill\"));
							circ[k].setAttributeNS(null, \"fill\", \"LightGray\");

						}
					}
				}
			curr_gene = tmp_curr_gene;

		}
		function showGenomesForNode(gene_id)
		{
			if (nodeIsClicked != 0)
			{
				var oldID = nodeIsClicked
				nodeIsClicked = 0;
				restoreNode(oldID);


			}
			if (gene_id != null && gene_id in node2genome)
			{
				var insertTd = document.getElementById(\"gene_list\");
				insertTd.innerHTML = \"\";
				var genomeNewList = node2genome[gene_id].split(\",\");
				var genomeNewStr = \"\";
				var curr_x = 0;
				for (var i =0; i < genomeNewList.length; i++)
				{
					if (curr_x + genomeNewList[i].length > 25)
					{
						genomeNewStr = genomeNewStr + \"<br />\";
						curr_x = 0;
					}
					genomeNewStr = genomeNewStr + genomeNewList[i] + \",\";
					curr_x = curr_x + genomeNewList[i].length + 1;
				}
				var circ =  document.getElementsByClassName(gene_id+\"2\");
				if (circ[0] != null)
				{
					circ[0].setAttributeNS(null, \"oldFill\", circ[0].getAttributeNS(null, \"fill\"));
					circ[0].setAttributeNS(null, \"fill\", \"green\");
					var oldR = circ[0].getAttributeNS(null, \"r\");
					circ[0].setAttributeNS(null, \"r\", oldR * 2);

				}
				insertTd.innerHTML = genomeNewStr;
				var insertPD = document.getElementById(\"percDiv\");
				insertPD.innerHTML = \"\";
				var str = \"\";

					var tot = null;
					var cnt = null;
					var circ = document.getElementsByClassName(gene_id + \"2\");
					if (circ != null)
					{
						for (var k =0 ;k < circ.length; k++)
						{
							if (circ[k].hasAttribute(\"totCnt\"))
							{
								tot = circ[k].getAttribute(\"totCnt\");
								cnt = circ[k].getAttribute(\"genomeCnt\");
							}
						}
					}
					if (tot != null)
					{
						str = str + \"Count:\" + cnt + \" (\" + parseInt(10000*cnt/tot)/100 +\")%\";
					}
					for (var i in metadata_is_checked)
					{
						if (metadata_is_checked[i])
						{
							var metaCnt = {};
							for (var j =0; j < genomeNewList.length; j++)
							{
								if (metaCnt !== undefined && metadata[genomeNewList[j]] !== undefined && metadata[genomeNewList[j]][i] in metaCnt)
								{
									metaCnt[metadata[genomeNewList[j]][i]] = metaCnt[metadata[genomeNewList[j]][i]] + 1;
								}
								else
								{
									if (metadata[genomeNewList[j]] !== undefined)
									{
										metaCnt[metadata[genomeNewList[j]][i]] = 1;
									}
								}
							}

							str = str + \"<br />\" + i + \": \";
							for (j in metaCnt)
							{
								str = str + j + \"(\" + metaCnt[j] + \"),\";
							}
							str.slice(0,-1);
						}
					}
				insertPD.innerHTML = str;

			}
		}

		//Restores a phylogeny node, the genome list and the percentage metadata to the default state
		function restoreNode(gene_id)
		{
			if (gene_id != null && gene_id in node2genome && nodeIsClicked== 0)
			{
				var insertTd = document.getElementById(\"gene_list\");
				var circ =  document.getElementsByClassName(gene_id+\"2\");
				if (circ[0] != null)
				{
					circ[0].setAttributeNS(null, \"fill\", circ[0].getAttributeNS(null, \"oldFill\"));
					var oldR = circ[0].getAttributeNS(null, \"r\");
					circ[0].setAttributeNS(null, \"r\", oldR / 2);

				}
				insertTd.innerHTML = \"\";
				var insertPd = document.getElementById(\"percDiv\");
				insertPd.innerHTML = \"\";
			}
		}

		//Makes an svg of a tree from the JSON file
		function draw_tree(level, type)
		{
			tree_level = level;
			curr_type = type; //This is linear and circular
			var meta_cnt = 0; //This counts the metadata variables turned on to get to appropiate height
			for (var i in metadata_is_checked)
			{
				if (metadata_is_checked[i] == true)
				{
					meta_cnt = meta_cnt + 1;
				}
			}

			//Make sure that there is a phylogeny loaded
			if (treeJSON != null)
			{

				//Reset the treeSVG
				var treeSVG = document.getElementById(\"tree_svg\");
				treeSVG.innerHTML=\"\";
				var height = treeSVG.getAttribute(\"height\");
				var width = treeSVG.getAttribute(\"width\");

				//Drawing a linear style tree
				if (type.toLowerCase() == \"linear\")
				{
					var y_dist = parseFloat(height); //Height is the y-distance
					var levs =  Object.keys(treeJSON).length +1 + meta_cnt; //Number of levels
					var x_dist = parseFloat(width)/(levs - parseInt(level)); //width over the number of levels is the x distance
					//i is the level of the trees
					for (var i in treeJSON)
					{
						//Only draw those levels greater than the level
						if (parseInt(i) >= level)
						{
							//Going through all branchs and nodes in the tree
							for (var j = 0; j < treeJSON[i].length; j++)
							{
								//this is vertical line
								if (treeJSON[i][j][6] ==0)
								{
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
									newSVG.setAttributeNS(null, \"class\",  treeJSON[i][j][7]+\"0\");
									newSVG.setAttributeNS(null, \"x1\", ((parseFloat(treeJSON[i][j][0])-level)*x_dist));
									newSVG.setAttributeNS(null, \"y1\", ((parseFloat(treeJSON[i][j][1])*y_dist)));
									newSVG.setAttributeNS(null, \"x2\", ((parseFloat(treeJSON[i][j][2])-level)*x_dist));
									newSVG.setAttributeNS(null, \"y2\", ((parseFloat(treeJSON[i][j][3]))*y_dist));
									newSVG.setAttributeNS(null, \"stroke\", \"black\");
									newSVG.setAttributeNS(null, \"stroke-width\", \"2\");
									treeSVG.appendChild(newSVG);
								}
								//this is horizontal lines
								if (treeJSON[i][j][6] == 1)
								{
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
									newSVG.setAttributeNS(null, \"class\",  treeJSON[i][j][7]+\"1\");
									newSVG.setAttributeNS(null, \"x1\", ((parseFloat(treeJSON[i][j][0])-level)*x_dist));
									newSVG.setAttributeNS(null, \"y1\", ((parseFloat(treeJSON[i][j][1])*y_dist)));
									newSVG.setAttributeNS(null, \"x2\", ((parseFloat(treeJSON[i][j][2])-level)*x_dist));
									newSVG.setAttributeNS(null, \"y2\", ((parseFloat(treeJSON[i][j][3]))*y_dist));
									newSVG.setAttributeNS(null, \"stroke\", \"black\");
									newSVG.setAttributeNS(null, \"stroke-width\", \"2\");
									treeSVG.appendChild(newSVG);
								}
								//this is the nodes
								if (treeJSON[i][j][6] == 2)
								{
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
									newSVG.setAttributeNS(null, \"class\",  treeJSON[i][j][7]+\"2\");
									newSVG.setAttributeNS(null, \"cx\", ((treeJSON[i][j][0]-level)*x_dist));
									newSVG.setAttributeNS(null, \"cy\", ((treeJSON[i][j][1])*y_dist));
									newSVG.setAttributeNS(null, \"r\", (y_dist/(treeJSON[i][j][5] * 3)));
									newSVG.setAttributeNS(null, \"cnt\",  treeJSON[i][j][4]);
									newSVG.setAttributeNS(null, \"fill\", \"black\");
									newSVG.setAttributeNS(null, \"onmouseover\", \"showGenomesForNode(\'\" + treeJSON[i][j][7] + \"\')\");
									newSVG.setAttributeNS(null, \"onmouseout\", \"restoreNode('\" + treeJSON[i][j][7] + \"\')\");
									newSVG.setAttributeNS(null, \"onclick\", \"keepNodeOn(\'\"+ treeJSON[i][j][7] + \"\')\");

									treeSVG.appendChild(newSVG);

								}
							}
						}
					}

				}
				//Drawing a circular style tree using the same
				if (type.toLowerCase() == \"circular\")
				{
					//xc is center node x-value, yc is the
					var xc = parseFloat(width) /2; var yc = parseFloat(height) /2;
					var r = xc;
					if (yc < r)
					{
						r = yc;
					}
					var levs =  Object.keys(treeJSON).length;
					var x_dist = r/(meta_cnt + levs - parseInt(level));
					var pi = Math.PI;
					for (var i in treeJSON)
					{
						if (parseInt(i) >= level)
						{
							for (var j = 0; j < treeJSON[i].length; j++)
							{
								//Getting the x & y locations for the \"polar\" values
								var x1 = Math.cos(parseFloat(treeJSON[i][j][1]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][0])) + xc;
								var y1 = Math.sin(parseFloat(treeJSON[i][j][1]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][0])) + yc;
								var x2 = Math.cos(parseFloat(treeJSON[i][j][3]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][2])) + xc;
								var y2 = Math.sin(parseFloat(treeJSON[i][j][3]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][2])) + yc;

								//Making the polar trees for each level
								//0 is the vertical branches
								//1 is the horizontal branches
								//2 is the node circle
								if (treeJSON[i][j][6] == 0)
								{
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
									newSVG.setAttributeNS(null, \"class\", treeJSON[i][j][7] + \"0\");
									newSVG.setAttributeNS(null, \"x1\", x1);
									newSVG.setAttributeNS(null, \"y1\", y1);
									newSVG.setAttributeNS(null, \"x2\", x2);
									newSVG.setAttributeNS(null, \"y2\", y2);
									newSVG.setAttributeNS(null, \"stroke\", \"black\");
									newSVG.setAttributeNS(null, \"stroke-width\", \"2\");

									treeSVG.appendChild(newSVG);
								}

								if (treeJSON[i][j][6] == 1)
								{
									var rad = x_dist * (levs - parseFloat(treeJSON[i][j][2]));
									var dir = 0;
									if (parseFloat(treeJSON[i][j][3]) > parseFloat(treeJSON[i][j][1]))
									{
										dir = \"1\";
									}
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"path\");
									newSVG.setAttributeNS(null, \"class\", treeJSON[i][j][7] + \"1\");
									newSVG.setAttributeNS(null, \"d\", \"M\" + x1+\" \" +y1+\" A \" + rad + \" \" + rad + \" 0 0 \" +dir+\" \" +x2+\" \" + y2);
									newSVG.setAttributeNS(null, \"stroke\", \"black\");
									newSVG.setAttributeNS(null, \"stroke-width\", \"2\");
									newSVG.setAttributeNS(null, \"fill\", \"none\");
									treeSVG.appendChild(newSVG);
								}

								if (treeJSON[i][j][6] == 2)
								{
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
									newSVG.setAttributeNS(null, \"class\", treeJSON[i][j][7] + \"2\");
									newSVG.setAttributeNS(null, \"cnt\",  treeJSON[i][j][4]);
									newSVG.setAttributeNS(null, \"cx\", x1);
									newSVG.setAttributeNS(null, \"cy\", y1);
									newSVG.setAttributeNS(null, \"r\", (x_dist/3));
									newSVG.setAttributeNS(null, \"fill\", \"black\");
									newSVG.setAttributeNS(null, \"stroke\", \"0\");
								//Node can have a mouse over and click functionality
									newSVG.setAttributeNS(null, \"onmouseover\", \"showGenomesForNode(\'\" + treeJSON[i][j][7] + \"\')\");
									newSVG.setAttributeNS(null, \"onmouseout\", \"restoreNode('\" + treeJSON[i][j][7] + \"\')\");
									newSVG.setAttributeNS(null, \"onclick\", \"keepNodeOn(\'\"+ treeJSON[i][j][7] + \"\')\");
									treeSVG.appendChild(newSVG);
								}
							}
						}
						var curr_x = meta_cnt + parseInt(level)-1;
						//Adding in the metadata values at the end of levels as they are clicked
						for (var k in metadata_head)
						{
							if (metadata_is_checked[k] == true)
							{
								curr_x = curr_x - 1;

								//Goes through the top level
								for (var j = 0; j < treeJSON[parseInt(level)].length; j++)
								{
									//Goes through only the vertical branches
									if (treeJSON[parseInt(level)][j][6] == 0)
									{
										var cx = Math.cos(parseFloat(treeJSON[parseInt(level)][j][3]) * 2 * pi) * x_dist * (levs - curr_x) + xc;
										var cy = Math.sin(parseFloat(treeJSON[parseInt(level)][j][3]) * 2 * pi) * x_dist * (levs - curr_x) + yc;
										var r = x_dist / 3;
										if (treeJSON[parseInt(level)][j][7] in node2meta)
										{
											var node_counts = node2meta[treeJSON[parseInt(level)][j][7]][k];
											var tot = treeJSON[parseInt(level)][j][4];
											var curr_p = 0;

											//Drawing the metadata circle/pie chart
											for (var i in node_counts)
											{
												var largeArc = 0;
												var perc = parseFloat(node_counts[i]);
												if ((perc / tot) > 0.49)
												{
													largeArc = 1; //To allow for arcs > 180 degrees
												}
												if (perc != tot)
												{
													var x1 = Math.cos(curr_p) * r + cx;
													var y1 = Math.sin(curr_p) * r + cy;
													var x2 = Math.cos((curr_p + perc / tot) * 2 * Math.PI) * r + cx;
													var y2 = Math.sin((curr_p + perc / tot) * 2 * Math.PI) * r + cy;
													curr_p = curr_p + perc / tot;
													var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"path\");
													newSVG.setAttributeNS(null, \"id\", j + \"3\");
													newSVG.setAttributeNS(null, \"d\", \"M\" + cx+\" \" +cy+\" L \" + x1+\" \" +y1+\" A \" + r + \" \" + r + \" 0 \"+ largeArc+ \" 1 \" +x2+\" \" + y2 + \" L \" + cx + \" \" + cy);
													newSVG.setAttributeNS(null, \"stroke\", \"red\");
													newSVG.setAttributeNS(null, \"stroke-width\", \"0\");
													newSVG.setAttributeNS(null, \"fill\", metadata_head[k][i]);

													treeSVG.appendChild(newSVG);
												}
												else
												{
													var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
													newSVG.setAttributeNS(null, \"id\", j + \"3\");
													newSVG.setAttributeNS(null, \"cx\", cx);
													newSVG.setAttributeNS(null, \"cy\", cy);
													newSVG.setAttributeNS(null, \"r\", r);
													newSVG.setAttributeNS(null, \"fill\", metadata_head[k][i]);

													treeSVG.appendChild(newSVG);
												}
											}
										}
										else
										{
											var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
											newSVG.setAttributeNS(null, \"id\", j + \"3\");
											newSVG.setAttributeNS(null, \"cx\", cx);
											newSVG.setAttributeNS(null, \"cy\", cy);
											newSVG.setAttributeNS(null, \"r\", r);
											newSVG.setAttributeNS(null, \"fill\", \"black\");

											treeSVG.appendChild(newSVG);

										}
									}
								}
							}
						}
					}
				}
				if (curr_gene != null)
				{
					add_to_tree(curr_gene);
				}
				if (nodeIsClicked != 0)
				{
					var oldNode = nodeIsClicked;
					showGenomesForNode(nodeIsClicked);
					nodeIsClicked = oldNode;
				}
			}
		}

		//Makes the fGI image-> used in resizes, etc
		function makeImage()
		{
			//Getting all the divisions and SVGs
			var divLeft = document.getElementById(\"leftDiv\");
			var divRight = document.getElementById(\"rightDiv\");
			var divLeg = document.getElementById(\"legendDiv\");
			var svgLeg = document.getElementById(\"legendSVG\");
			var svgLeg = document.getElementById(\"legendSVG\");

			var divMain = document.getElementById(\"mainDiv\");
			var divSave = document.getElementById(\"saveDiv\");
			var svgSave = document.getElementById(\"diskSVG\");
			var divSort = document.getElementById(\"sortDiv\");
			var svgMain = document.getElementById(\"mainSVG\");
			var svgBottom = document.getElementById(\"bottomRow\");
			var legRect = document.getElementById(\"legRect\");

			//Height of the footer to help set the height of the main image
			var footTextHeight = parseFloat(svgBottom.style.height)+parseFloat(svgBottom.style.top);
			svgMain.style.height = footTextHeight + \"px\"

			var headTextWidth = parseFloat(document.getElementById(\"headerText\").getBBox().width);
			var footTextWidth = parseFloat(document.getElementById(\"footerText\").getBBox().width);

			//Getting all the heights and widths of the different divisions and SVGs
			//These are the ones stored by the program, not calculated
			var totWidth = parseFloat(divMain.getAttribute(\"tot_width\"));
			var totHeight = parseFloat(divMain.getAttribute(\"tot_height\"));
			var leftWidth = parseFloat(divLeft.getAttribute(\"wd\"));
			var legWidth = parseFloat(divLeg.getAttribute(\"wd\"));
			var legHeight = parseFloat(divLeg.getAttribute(\"ht\"));
			var rightWidth = parseFloat(divRight.getAttribute(\"wd\"));
			var saveHeight = parseFloat(divSave.getAttribute(\"ht\"));
			var saveWidth = parseFloat(divSave.getAttribute(\"wd\"));

			//set the legend hiehg
			if (legHeight < saveHeight)
			{
				legHeight = saveHeight;
			}
			if (leftWidth < saveWidth)
			{
				leftWidth = saveWidth;
			}




			legRect.setAttribute(\"width\", totWidth-rightWidth-leftWidth);

			//(re-)Making the legend
			var curX = parseFloat(legRect.getAttribute(\"x\"))+10;
			var curY = 10;
			var maxY = 0;
			var maxX = 0;
			for (var i = 0; i < divLeg.getAttribute(\"num\"); i++)
			{
				var tmpText = document.getElementById(\"TextID\"+i);
				//bb is bounding box and gives the width and height of each string
				var bbWidth = tmpText.getBBox().width;
				var bbHeight = tmpText.getBBox().height;
				if (bbHeight > maxY)
				{
					maxY = bbHeight;
				}
				//Either extends the legend width or extends the row
				if (curX + bbWidth > legRect.getAttribute(\"width\"))
				{

					if (maxX < curX)
					{
						maxX = curX;
					}
					curX = parseFloat(legRect.getAttribute(\"x\"))+10;
					curY = curY + maxY + 1;
					MaxY = 0;
				}
				tmpText.setAttribute(\"x\", curX);
				tmpText.setAttribute(\"y\", curY);
				curX = curX + bbWidth+10;

			}
			//Resetting maxX as needed
			if (curX > maxX)
			{
				maxX = curX;
			}

			//Resetting the Legend Box
			var tmpText = document.getElementById(\"LegText\");
			var bbWidth = tmpText.getBBox().width;
			var bbHeight = tmpText.getBBox().height;
			tmpText.setAttribute(\"y\", curY + maxY + bbHeight);
			legHeight = bbHeight + curY +maxY + 10;
			legWidth = maxX + 10;
			legRect.setAttribute(\"height\", legHeight);
			legRect.setAttribute(\"width\", legWidth);
			tmpText.setAttribute(\"x\", bbWidth/2 + (legWidth - bbWidth)/2);

			//Resetting division heights and widths
			if (divSort.getBoundingClientRect().height > divSort.style.height)
			{
				divSort.style.height = divSort.getBoundingClientRect().height;
			}
			if (divSave.getBoundingClientRect().height > divSave.style.height)
			{
				divSave.style.height = divSave.getBoundingClientRect().height;
			}

			if (divSort.style.height > legHeight)
			{
				legHeight = divSort.style.height;
			}
			if (divSave.style.height > legHeight)
			{
				legHeight = divSave.style.height;
			}
			divLeg.style.height = legHeight;
			divLeg.style.width = legWidth;

			//Getting the total window size
			var windWidth = window.innerWidth;
			var windHeight = window.innerHeight;


			//Re-calculating totHeight and totWidth where neaded
			if (parseFloat(svgMain.style.width) < headTextWidth)
			{
				svgMain.style.width = headTextWidth;
			}
			if (parseFloat(svgMain.style.width) < footTextWidth)
			{
				svgMain.style.width = footTextWidth;
			}
			if (parseFloat(svgMain.style.height) < parseFloat(document.getElementById(\"bottomRow\").style.height) + parseFloat(document.getElementById(\"bottomRow\").style.top))
			{
				svgMain.style.height = parseFloat(document.getElementById(\"bottomRow\").style.height) + parseFloat(document.getElementById(\"bottomRow\").style.top);
			}
			if (totHeight < parseFloat(svgMain.style.height) + parseFloat(divLeg.style.height))
			{
				totHeight = parseFloat(svgMain.style.height) + parseFloat(divLeg.style.height)
			}

			if (parseFloat(divMain.scrollHeight) > parseFloat(divMain.clientHeight))
			{
				totWidth = totWidth + 10;
			}

			if (totWidth < 10+leftWidth + rightWidth +parseFloat(svgMain.style.width))
			{
				totWidth = leftWidth + rightWidth +parseFloat(svgMain.style.width)+10;
			}
			if (totWidth > windWidth)
			{
				totWidth = windWidth;
			}

			if (totHeight > windHeight)
			{
				totHeight = windHeight;
				rightWidth = rightWidth+40;
				totWidth = totWidth+80;
				if (totWidth > windWidth-40)
				{
					totWidth = windWidth-40;
				}

			}

			//Using total window size to reset all the divisions/svgs
			divLeft.style.position = \"absolute\";
			divLeft.style.overflowY = \"hidden\";

			divLeft.style.top = 0;
			divLeft.style.left = 0;
			divLeft.style.width = leftWidth;
			divLeft.style.height = (totHeight - legHeight);

			//divSave is bottom right
			divSave.style.top = totHeight-legHeight;
			divSave.style.left = 0;
			divSave.style.width = leftWidth;
			divSave.style.height = legHeight;
			divSave.style.position = \"absolute\";

			svgSave.style.top = 0;
			svgSave.style.left = 0;
			svgSave.style.width = leftWidth;
			svgSave.style.position = \"absolute\";



			divMain.style.position = \"absolute\";
			divMain.style.overflowY = \"auto\";
			divMain.style.overflowX = \"auto\";
			divMain.style.top = 0;
			divMain.style.left = leftWidth;
			divMain.style.width = (totWidth-rightWidth-leftWidth-10);
			divMain.style.height = (totHeight - legHeight);

			divRight.style.position = \"absolute\";
			divRight.style.overflowY = \"hidden\";
			divRight.style.scrollWidth =0;
			divRight.style.top = 0;
			divRight.style.left = (totWidth-rightWidth);
			divRight.style.width = rightWidth;
			divRight.style.height = totHeight - legHeight;

			divSort.style.position = \"absolute\";
			divSort.style.top = totHeight-legHeight;
			divSort.style.left = (totWidth-rightWidth);
			divSort.style.width = rightWidth;
			divSort.style.height =legHeight;

			divLeg.style.overflowX = \"auto\";
			divLeg.style.position = \"absolute\";
			divLeg.style.top = totHeight-legHeight;
			divLeg.style.left= leftWidth;
			divLeg.style.width = totWidth-rightWidth-leftWidth;
			divLeg.style.height = legHeight;
			svgLeg.style.height = legHeight-5;
			svgLeg.style.top = 0;
			svgLeg.style.left= (totWidth-rightWidth-leftWidth-maxX-10)/2;
			svgLeg.style.width = maxX+10;

			if (parseFloat(document.getElementById(\"topRow\").style.width) !=  parseFloat(svgMain.style.width))
			{
				document.getElementById(\"topRow\").style.width = parseFloat(svgMain.style.width);
			}

			if (parseFloat(document.getElementById(\"bottomRow\").style.width) !=  parseFloat(svgMain.style.width))
			{
				document.getElementById(\"bottomRow\").style.width = parseFloat(svgMain.style.width);
			}
			svgLeg.style.position = \"absolute\";
			document.getElementById(\"footerText\").setAttribute(\"x\", parseFloat(divMain.style.width)/2);
			document.getElementById(\"headerText\").setAttribute(\"x\", parseFloat(divMain.style.width)/2);

		}

		//Function that shows/hides all the singletons
		function showSingleton()
		{
			var svgMain = document.getElementById(\"mainSVG\");
			var topRow = document.getElementById(\"topRow\");
			var bottomRow = document.getElementById(\"bottomRow\");
			var rowNum = parseInt(svgMain.getAttribute(\"numrow\"));
			var divSort = document.getElementById(\"sortDiv\");

			//type is whether the singletons are on (1) or off (0)
			var type = divSort.getAttribute(\"singleton\");
			var names = JSON.parse(namesJSON);
			var cur_y = 0;

			//going through all the fGIs and checking for the fGIs that have only one genome
			//these are the singletons
			for (var i = 0; i < rowNum; i++)
			{
				var curRow = document.getElementById(\"mainRow\"+i);
				var rowID = curRow.getAttribute(\"rowid\");
				var curLeft = document.getElementById(\"rowLeft\"+i);
				var curRight = document.getElementById(\"All\"+names[i]);

				//getting the genome counts
				var cnt = parseInt(curRow.getAttribute(\"cnt\"));
				if (cnt <= 1)
				{
					//Turning the singletons off by hiding them
					if (type == \"1\")
					{
						curRow.setAttribute(\"visibility\", \"hidden\");
						curRow.setAttribute(\"y\", cur_y);
						curLeft.setAttribute(\"visibility\", \"hidden\");
						curLeft.style.top = cur_y;
						curRight.style.visibility\ = \"hidden\";
						curRight.style.top = cur_y;
					}
					else
					{
						//Turning the singletons on by making them visible
						curRow.setAttribute(\"visibility\", \"visible\");
						curRow.setAttribute(\"y\", cur_y);
						curLeft.setAttribute(\"visibility\", \"visible\");
						curLeft.style.top = cur_y;
						curRight.style.visibility=  \"visible\";
						curRight.style.top = cur_y;
						cur_y = cur_y + parseFloat(curRow.getAttribute(\"height\"));

					}
				}
				else
				{
					//Restting y-position as needed
					curRow.setAttribute(\"y\", cur_y);
					curRight.style.top =  cur_y;
					curLeft.style.top = cur_y;

					cur_y = cur_y + parseFloat(curRight.getAttribute(\"height\"));
				}
			}
			bottomRow.style.top = cur_y + parseFloat(topRow.style.height);
			svgMain.style.height = cur_y;

			//Changes singleton hide and show button as needed
			if (type == \"1\")
			{
				divSort.setAttribute(\"singleton\", \"0\");
				var singText = document.getElementById(\"singletonButton\");
				singText.innerHTML = \"Add Singletons\";
			}
			else
			{
				divSort.setAttribute(\"singleton\", \"1\");
				var singText = document.getElementById(\"singletonButton\");
				singText.innerHTML = \"Remove Singletons\";
			}
		}

		//Change the save type (PNG or SVG) button
		function changeOut(type)
		{
			var typeNode = document.getElementById(\"rect\" + type);
			var oldNode = document.getElementById(\"rect\" + saveType);
			if (typeNode != null)
			{
				oldNode.setAttribute(\"fill\", \"#bbbbbb\");
				typeNode.setAttribute(\"fill\", \"#bb00bb\");
				saveType = type;
			}
		}

		//Saves the multiple fasta of the selected fGIs to a text file
		function saveMultiFasta()
		{
			var allFasta = JSON.parse(seqsJSON);
			var html = \"<p>\";
			var txt = \"\";
			for (var i in allFasta)
			{

					var outVar = allFasta[i];
					outVar1 = outVar.replace(/--ret--/g, \"\\n\");
					outVar2 = outVar.replace(/--ret--/g, \"<p>\");

					var fastaHeader1 = \">\";
					fastaHeader1 = fastaHeader1 + i + \"\\n\";

					var fastaHeader2 = \"&gt;\";
					fastaHeader2 = fastaHeader2 + i + \"<p>\";

					txt = txt + fastaHeader1 + outVar1;
					html = html + fastaHeader2 + outVar2;
			}
			html = html + \"</p>\";

			var nw = document.createElement(\"a\");
			var inTxt = new Blob([txt], {type: \'text\'});
			var txtfile = URL.createObjectURL(inTxt);
			var location = window.location.href;
			var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);
			nw.setAttribute(\"download\", currentID +\" .all.fasta\");
			nw.setAttribute(\"href\", txtfile);
			nw.setAttribute(\"target\", \"_blank\");
			alert(\"Saved TXT to \" + currentID +\".all.fasta\");
			document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
		}

		//Makes a single gene page with multiple tabs
		//used in both fGI and Core genes
		function writeFasta(id)
		{
			var fasta = JSON.parse(fastaJSON);
			var type = JSON.parse(outJSON);
			var terms = JSON.parse(termJSON);

			//If the gene passed along ID, write the html
			if (id in fasta)
			{
				var html = \"<html><script type=\'text/javascript\' src=\'../json/\"+id+\".allFasta.json\'></script><script type=\'text/javascript\' src=\'https\://s3.eu-central-1.amazonaws.com/cdn.bio.sh/msa/latest/msa.min.gz.js\'></script>\";
				if (treeJSON != null)
				{
					html = html + \"<script type=\'text/javascript\' src=\'../json/tree.json\'></script>\";
				}
				html = html + \"<script type=\'text/javascript\' src=\'../scripts/geneFile.functions.js\'></script><body style=\'font-family:Courier New\'>\";
				//Adding tab radio buttons
				html = html + \"<div id=\'select\'><form><input type=\'radio\' name=\'display\' value=\'summary\' checked onclick=\'changeDiv(event)\'>Summary</input>\";
				html = html + \"<input type=\'radio\' name=\'display\' value=\'centFasta\' onclick=\'changeDiv(event)\'>Centroid Fasta</input>\";
				html = html + \"<input type=\'radio\' name=\'display\' value=\'multFastaView\' onclick=\'changeDiv(event)\'>Multi-Fasta Viewer</input>\";
				html = html + \"<input type=\'radio\' name=\'display\' value=\'allFasta\' onclick=\'changeDiv(event)\'>Download Multi-Fasta File</input>\";
				if (treeJSON != null)
				{
					html = html + \"<input type=\'radio\' name=\'display\' value=\'tree\' onclick=\'changeDiv(event)\'>View Phylogeny</input>\";

				}
				//Adding in the summary information in a table
				html = html + \"</form></div><div id=\'summary\' style=\'visibility:visible;\'><table >\";
				for (var j in type)
				{
					var i = type[j];
					html = html + \"<tr><td style=\'border:1px solid black;\' id=\'ID\" + i + \"\'>\"+i+\"</td>\";
					var outVar = fasta[id][i];
					if (i == \"Sequence\")
					{
						outVar = outVar.replace(/--ret--/g, \"<br>\");
					}
					if (i == \"Associated Terms\")
					{
						var tmpVar = outVar.split(\";\");
						outVar = \"\";
						for (j in tmpVar)
						{
							if (tmpVar[j] in terms)
							{
								outVar = outVar + tmpVar[j] +\"(\"+ terms[tmpVar[j]] +\")<br>\";
							}
						}
					}
					html = html + \"<td style=\'border:1px solid black;\' id=\'\" + i + \"Value\'>\"+outVar+\"</td></tr>\";
				}
				html = html + \"</table></div><div id=\'centFasta\' style=\'visibility:hidden;\'>\";
				//adding in the fasta pages
				var outVar = fasta[id][\"Sequence\"];
				outVar = outVar.replace(/--ret--/g, \"<br>\");
				var fastaHeader = \"&gt;\";
				fastaHeader = fastaHeader + \"Centroid_\" + id + \"<br>\";
				html = html + fastaHeader + outVar + \"</p></div><div id=\'multFastaView\' style=\'visibility:hidden;\'><div id=\'multFastaViewInner\'></div>Viewer is supported by <a href=\'msa.biojs.net\'>MSA Viewer</a></div><div id=\'allFasta\' style=\'visibility:hidden;\'>\";
				//Making tree page if nessicary
				if (treeJSON != null)
				{
					html = html + \"</div><div id='tree' style=\'visibility:hidden;\'>\";

					html = html + \"<table ><tr height=\\\"\" + (tree_size * 0.2) + \"\\\"><td style=\\\"vertical-align:top;\\\" rowspan=\\\"3\\\">\"

					html = html + \"<div><svg version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" id=\\\"tree_svg\\\" height=\\\"\"+ tree_size+\"\\\" width=\\\"\"+ tree_size+\"\\\"></svg>\";

					html = html + \"</div></td><td style=\\\"vertical-align:top;\\\"><div><h3>Tree Functions</h3>Outer Ring Tree Level:<button onclick=\\\"changeTree(\'-\')\\\">-</button><b id = \\\"LevelCount\\\">4</b><button onclick=\\\"changeTree(\'+\')\\\">+</button>\";

					html = html + \"<form>Phylogeny Style:<input type=\\\"radio\\\" name=\\\"type\\\" value=\\\"Circular\\\" checked=\\\"true\\\" onclick=\\\"changeType(event)\\\"/>Circular<input type=\\\"radio\\\" name=\\\"type\\\" value=\\\"Linear\\\" onclick=\\\"changeType(event)\\\"/>Linear</form>\";

					html = html + \"<svg version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" id=\\\"save_svg\\\" height=\\\"\"+ (tree_size *0.1)+\"\\\" width=\\\"\"+ (tree_size*0.4)+\"\\\"><text x=\\\"2\\\" y=\\\"14\\\">Node Colors:</text><circle cx=\\\"155\\\" r=\\\"4\\\" cy=\\\"8\\\" fill=\\\"black\\\"/><text x=\\\"162\\\" y=\\\"14\\\">Present</text><circle cx=\\\"245\\\" r=\\\"4\\\" cy=\\\"8\\\" fill=\\\"LightGray\\\"/><text x=\\\"252\\\" y=\\\"14\\\">Missing</text>$disk_image_svg $disk_image_png</svg></div></td></tr><tr height=\\\"\"+ (tree_size *0.4)+\"\\\" >\";
					html = html + \"<td style=\\\"vertical-align:top;\\\"><div id=\\\"gene_list\\\" style=\\\"height:\"+ (tree_size*0.4)+\"; overflow-y: auto;\\\"></div></td></tr><tr height=\\\"\"+ (tree_size*0.15)+\"\\\"><td style=\\\"vertical-align:bottom;\\\">\";
					html = html + \"<div id=\\\"percentage_image\\\"></div></td></tr><tr height=\\\"\"+ (tree_size*0.25)+\"\\\">\";
					html = html + \"<td><table><tr><td id=\\\"SelectMetatype\\\" width=\\\"\"+ (tree_size*0.4)+\"\\\" style=\\\"vertical-align:top;\\\"><h3>Metadata Type</h3>\";
					for (var a in metadata_head)
					{
						html = html + \"<input type=\\\"checkbox\\\" id=\\\"select\"+a+\"\\\" value=\\\"\" + a+ \"\\\" onchange=\\\"change_metadata_status(event)\\\" />\"+a+\"<p>\";
					}
					html = html +\"</td><td id = \\\"Legend\\\" style=\\\"vertical-align:top;\\\"><h3>Legend</h3><svg id=\\\"LegendSVG\\\" version=\\\"1.2\\\" overflow=\\\"visible\\\" baseProfile=\\\"tiny\\\" xmlns=\\\"http://www.w3.org/2000/svg\\\" height=\\\"\"+ (tree_size*0.2)+\"\\\" width=\\\"\"+ (tree_size*0.6)+\"\\\"></svg></td></tr></table><td style=\\\"vertical-align:top;\\\"><h3>Percentages</h3><div id =\\\"percDiv\\\"></div></td></tr></table>\";
				}
				html = html +\"</div></body></html>\";
				//Writing html to a new page
				var newWind = window.open();
				newWind.document.write(html);
				newWind.document.close();

			}
			else
			{
				alert(\"No Sequence information for \" + id);
			}
		}

	//Save the Tree as an SVG or PNG: used for core region page
	function saveSVGtree(fileType, id)
	{
		var event = new Event(\'build\');
		var nw = document.createElement(\"a\");
		var svg_all = document.getElementById(\"tree_svg\");
		var add = null;
		var FGRmodel = document.getElementById(\"selectFGR\");
		if (FGRmodel)
		{
			console.log(\"Working\");
			var rect = parseFloat(svg_all.getAttribute(\"height\"));
			add =  \"<g id = \\\"selectFGR\\\" transform=\\\"translate(0,\" + rect+ \")\\\">\"+ FGRmodel.innerHTML + \"</svg>\";
		
		}
		var svgStr = \"<svg xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" height=\\\"$SVGHEIGHT\\\" width=\\\"$SVGWIDTH\\\">\"+document.getElementById(\"LegendSVG\").innerHTML+ document.getElementById(\"tree_svg\").innerHTML + \"</svg>\";
		if (add != null)
		{
			svgStr = \"<svg xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" height=\\\"$SVGHEIGHT\\\" width=\\\"$SVGWIDTH\\\">\"+document.getElementById(\"LegendSVG\").innerHTML+ document.getElementById(\"tree_svg\").innerHTML + add+\"</svg>\";
		}
		var inSvg = new Blob([svgStr], {type: \'image/svg+xml\'});
		var svgfile = URL.createObjectURL(inSvg);
		var location = window.location.href;
		var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);
		if (fileType==\"png\")
		{
			var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"svg\");
			newSVG.setAttribute(\"xmlns:xlink\",\"http://www.w3.org/1999/xlink\");
			newSVG.setAttribute(\"height\", parseFloat(svg_all.getAttribute(\"width\")) * 2);
			newSVG.setAttribute(\"width\",2*parseFloat(svg_all.getAttribute(\"height\")) + 200);
			newSVG.innerHTML = svgStr;
			var svgStr = new XMLSerializer().serializeToString(newSVG);

			var canvas = document.createElement(\"canvas\");
			canvas.width = parseFloat(svg_all.getAttribute(\"width\")) * $DPI;
			canvas.height = (parseFloat(svg_all.getAttribute(\"height\")) + 200) * $DPI;
			var ctx = canvas.getContext(\"2d\");
			ctx.scale($DPI, $DPI);

			var img = new Image();
			img.src = \"data:image/svg+xml;base64,\" + window.btoa(svgStr);
			img.onload = function()
			{
				ctx.drawImage(img, 0, 0);

				nw.setAttribute(\"download\",id + \"tree.$filenamePNG\");
				nw.setAttribute(\"href\", canvas.toDataURL(\"image/png\"));
				nw.setAttribute(\"target\", \"_blank\");

				alert(\"Saved PNG to \" + id + \"tree.$filenamePNG\");
				document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
			
			}
		}
		else
		{
			nw.setAttribute(\"download\", id + \".tree.$filenameSVG\");
			nw.setAttribute(\"href\", svgfile);
			nw.setAttribute(\"target\", \"_blank\");
			alert(\"Saved SVG to \"  + id + \".tree.$filenameSVG\");
			document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
		}
	}


//This is used to seed the single gene summary page
var outJSON=\'[\"Cluster\", \"Name\", \"Region\", \"# of Genomes\",\"Functional IDs\",\"Associated Terms\", \"Mean\", \"Standard Deviation\", \"Minimum Length\", \"Maximum Length\", \"Sequence\", \"Genomes\"]\';
";

}

sub write_geneFile_javascript {

    open( FILE, ">", $out_dir . "/$output_id/scripts/geneFile.functions.js" );
    print FILE "
var onDiv = \"summary\"; //The tab shown on loading the single gene page
var allFasta = JSON.parse(allFastaJSON); //array with sequences
var html = \"<p>\"; //
var txt = \"\"; //Text to store the fasta string
var isMOn = 0; //flag on whether the multi-fasta viewer is on or off
var isTreeOn = 0; //whether the tree is showns or not
var tree_level; //level that is shown on edge. 1 are the leaves
var treeJSON = null; //whether the gene is
var genome2node; //
var node2genome = []; //
var metadata_head; //
var metadata; //array that maps metadata values to
var node2meta; //array that maps each node to the metadata values
var nodeIsClicked = 0; //which tree node is clicked
var meta_color = []; //array with the metadata color
var metadata_is_checked = []; //whether the different metadata is checked
var all_tree_nodes = []; //all the nodes in the tree
var geneList; //List of the genomes
var curr_genome = null;
var curr_gene = null;
var curr_type = \"circular\";

if (typeof tree_format !== \'undefined\')
{
	treeJSON = JSON.parse(tree_format);
	for (var i in treeJSON)
	{
		for (var j in treeJSON[i])
		{
			if (treeJSON[i][j][6] == \"2\")
			{
				all_tree_nodes[treeJSON[i][j][7]] = [];
				all_tree_nodes[treeJSON[i][j][7]][0] = i;
				all_tree_nodes[treeJSON[i][j][7]][1] = j;
			}
		}
	}
	genome2node = JSON.parse(gene2node);
	if (typeof meta_data != \'undefined\')
	{
		metadata_head = JSON.parse(meta_head_var);
		metadata = JSON.parse(meta_data);
		node2meta = JSON.parse(meta_node_data);
		for (var i in metadata_head)
		{
			metadata_is_checked[metadata_head[i]] = false;
		}
	}
	for (var i in genome2node)
	{
		var nodes = genome2node[i].split(\";\");
		for (var j = 0; j < nodes.length; j++)
		{
			if (!nodes[j] in node2genome || typeof node2genome[nodes[j]] == \'undefined\')
			{
				node2genome[nodes[j]] = i;
			}
			else
			{
				node2genome[nodes[j]] = node2genome[nodes[j]] + \",\" + i;
			}
		}
	}
}
for (var i =0; i < allFasta.length; i++)
{
	var outVar = allFasta[i][\"seq\"];
	outVar1 = outVar.replace(/--ret--/g, \"\\n\");
	outVar2 = outVar.replace(/--ret--/g, \"<p>\");

	var fastaHeader1 = \">\";
	fastaHeader1 = fastaHeader1 + allFasta[i][\"id\"] + \"\\n\";

	var fastaHeader2 = \"&gt;\";
	fastaHeader2 = fastaHeader2 + allFasta[i][\"id\"] + \"<p>\";

	txt = txt + fastaHeader1 + outVar1;
	html = html + fastaHeader2 + outVar2;
}
html = html + \"</p>\";


//Changing the Selected Tab in the
function changeDiv(evt)
{
	var target = evt.target;

	if (target.value != onDiv)
	{
		//saving the fasta file is selected
		if (target.value == \"allFasta\")
		{
			var nw = document.createElement(\"a\");
			var inTxt = new Blob([txt], {type: \'text\'});
			var txtfile = URL.createObjectURL(inTxt);
			var location = window.location.href;
			var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);
			nw.setAttribute(\"download\", currentID +\" .all.fasta\");
			nw.setAttribute(\"href\", txtfile);
			nw.setAttribute(\"target\", \"_blank\");
			alert(\"Saved TXT to \" + currentID +\".all.fasta\");
			document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
		}

		//If they change the tab to the multi-fasta viewer, load it if it is for the first time
		if (target.value == \"multFastaView\" && isMOn == 0)
		{
			var seqs =  msa.io.fasta.parse(txt);

			var m = msa({
				el: document.getElementById(\"multFastaViewInner\"),
				seqs: seqs,
				vis: {
					conserv: false,
					overviewbox: false
				},
				colorscheme: {scheme: \"nucleotide\"},
				menu: \"small\",
				bootstrapMenu: true
			});

			m.render(); //makes the viewer
			isMOn = 1;
		}
		if (target.value == \"tree\" && isTreeOn == 0)
		{
			draw_tree(1, curr_type);
			var genomes = document.getElementById(\"GenomesValue\").innerHTML;
			var new_list = genomes.replace(/<br>/g, \"\");
			 show_nodes_on_tree(new_list);
			 isTreeOn = 1;
		}
		var old = document.getElementById(onDiv);
		old.style.visibility = \"hidden\";
		onDiv = target.value;
		var nw = document.getElementById(onDiv);
		nw.style.visibility = \"visible\";
		var par = nw.parentNode;
		par.insertBefore(nw, old);
	}
}

//Draw the phylogeny on the tree page
function draw_tree(level, type)
{
	tree_level = level;
	curr_type = type;
	var meta_cnt = 0;
	for (var i in metadata_is_checked)
	{
		if (metadata_is_checked[i])
		{
					meta_cnt = meta_cnt + 1;
				}
			}
			if (treeJSON != null)
			{
				var treeSVG = document.getElementById(\"tree_svg\");
				treeSVG.innerHTML=\"\";
				var height = treeSVG.getAttribute(\"height\");
				var width = treeSVG.getAttribute(\"width\");
				if (type.toLowerCase() == \"linear\")
				{
					var y_dist = parseFloat(height);
					var levs =  Object.keys(treeJSON).length +1 + meta_cnt;
					var x_dist = parseFloat(width)/(levs - parseInt(level));
					for (var i in treeJSON)
					{
						if (parseInt(i) >= level)
						{
						for (var j = 0; j < treeJSON[i].length; j++)
						{
							if (treeJSON[i][j][6] ==0)
							{
								var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
								newSVG.setAttributeNS(null, \"class\",  treeJSON[i][j][7]+\"0\");
								newSVG.setAttributeNS(null, \"x1\", ((parseFloat(treeJSON[i][j][0])-level)*x_dist));
								newSVG.setAttributeNS(null, \"y1\", ((parseFloat(treeJSON[i][j][1])*y_dist)));
								newSVG.setAttributeNS(null, \"x2\", ((parseFloat(treeJSON[i][j][2])-level)*x_dist));
								newSVG.setAttributeNS(null, \"y2\", ((parseFloat(treeJSON[i][j][3]))*y_dist));
								newSVG.setAttributeNS(null, \"stroke\", \"black\");
								newSVG.setAttributeNS(null, \"stroke-width\", \"2\");
								treeSVG.appendChild(newSVG);
							}
							if (treeJSON[i][j][6] == 1)
							{
								var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
								newSVG.setAttributeNS(null, \"class\",  treeJSON[i][j][7]+\"1\");
								newSVG.setAttributeNS(null, \"x1\", ((parseFloat(treeJSON[i][j][0])-level)*x_dist));
								newSVG.setAttributeNS(null, \"y1\", ((parseFloat(treeJSON[i][j][1])*y_dist)));
								newSVG.setAttributeNS(null, \"x2\", ((parseFloat(treeJSON[i][j][2])-level)*x_dist));
								newSVG.setAttributeNS(null, \"y2\", ((parseFloat(treeJSON[i][j][3]))*y_dist));
								newSVG.setAttributeNS(null, \"stroke\", \"black\");
								newSVG.setAttributeNS(null, \"stroke-width\", \"2\");
								treeSVG.appendChild(newSVG);
							}
							if (treeJSON[i][j][6] == 2)
							{
								var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
								newSVG.setAttributeNS(null, \"class\",  treeJSON[i][j][7]+\"2\");
								newSVG.setAttributeNS(null, \"cx\", ((treeJSON[i][j][0]-level)*x_dist));
								newSVG.setAttributeNS(null, \"cy\", ((treeJSON[i][j][1])*y_dist));
								newSVG.setAttributeNS(null, \"r\", (y_dist/(treeJSON[i][j][5] * 3)));
								newSVG.setAttributeNS(null, \"cnt\",  treeJSON[i][j][4]);
								newSVG.setAttributeNS(null, \"fill\", \"black\");
								newSVG.setAttributeNS(null, \"onmouseover\", \"showGenomesForNode(\'\" + treeJSON[i][j][7] + \"\')\");
								newSVG.setAttributeNS(null, \"onmouseout\", \"restoreNode('\" + treeJSON[i][j][7] + \"\')\");
								newSVG.setAttributeNS(null, \"onclick\", \"keepNodeOn(\'\"+ treeJSON[i][j][7] + \"\')\");

								treeSVG.appendChild(newSVG);

							}
						}
						}
					}

				}
				if (type.toLowerCase() == \"circular\")
				{
					//Setting the polar tree variables:
					//xc = center x-axis location; yx = center y-axis location; r = radius size
					var xc = parseFloat(height) /2; var yc = parseFloat(width) /2;
					var r = xc;
					if (yc < r)
					{
						r = yc;
					}
					var levs =  Object.keys(treeJSON).length;
					var x_dist = r/(meta_cnt + levs - parseInt(level));
					var pi = Math.PI;
					//Tree format: made modularly with each node
					//3 types: 0 = vertical branches, 1 = horizontal, 2 = node
					//all are specifed in the treeJSON file that has been loaded
					for (var i in treeJSON)
					{
						//Only drawing the levels above the current level
						if (parseInt(i) >= level)
						{

							for (var j = 0; j < treeJSON[i].length; j++)
							{
								var x1 = Math.cos(parseFloat(treeJSON[i][j][1]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][0])) + xc;
								var y1 = Math.sin(parseFloat(treeJSON[i][j][1]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][0])) + yc;
								var x2 = Math.cos(parseFloat(treeJSON[i][j][3]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][2])) + xc;
								var y2 = Math.sin(parseFloat(treeJSON[i][j][3]) * 2 * pi) * x_dist * (levs - parseFloat(treeJSON[i][j][2])) + yc;
								//Making vertical branches- these are lines
								if (treeJSON[i][j][6] == 0)
								{
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"line\");
									newSVG.setAttributeNS(null, \"class\", treeJSON[i][j][7] + \"0\");
									newSVG.setAttributeNS(null, \"x1\", x1);
									newSVG.setAttributeNS(null, \"y1\", y1);
									newSVG.setAttributeNS(null, \"x2\", x2);
									newSVG.setAttributeNS(null, \"y2\", y2);
									newSVG.setAttributeNS(null, \"stroke\", \"black\");
									newSVG.setAttributeNS(null, \"stroke-width\", \"2\");

									treeSVG.appendChild(newSVG);
								}
								//Making horizontal branches- these are arcs
								if (treeJSON[i][j][6] == 1)
								{
									var rad = x_dist * (levs - parseFloat(treeJSON[i][j][2]));

									//dir is necessary to make sure that the arc is completed
									var dir = 0;
									if (parseFloat(treeJSON[i][j][3]) > parseFloat(treeJSON[i][j][1]))
									{
										dir = \"1\";
									}
									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"path\");
									newSVG.setAttributeNS(null, \"class\", treeJSON[i][j][7] + \"1\");
									newSVG.setAttributeNS(null, \"d\", \"M\" + x1+\" \" +y1+\" A \" + rad + \" \" + rad + \" 0 0 \" +dir+\" \" +x2+\" \" + y2);
									newSVG.setAttributeNS(null, \"stroke\", \"black\");
									newSVG.setAttributeNS(null, \"stroke-width\", \"2\");
									newSVG.setAttributeNS(null, \"fill\", \"none\");
									treeSVG.appendChild(newSVG);
								}
								//Making nodes- these are circles
								if (treeJSON[i][j][6] == 2)
								{

									var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
									newSVG.setAttributeNS(null, \"class\", treeJSON[i][j][7] + \"2\");
									newSVG.setAttributeNS(null, \"cx\", x1);
									newSVG.setAttributeNS(null, \"cy\", y1);
									newSVG.setAttributeNS(null, \"r\", (x_dist/3));
									newSVG.setAttributeNS(null, \"fill\", \"black\");
									newSVG.setAttributeNS(null, \"stroke\", \"0\");
									newSVG.setAttributeNS(null, \"cnt\",  treeJSON[i][j][4]);

									//Turning on the functionality of the nodes
									newSVG.setAttributeNS(null, \"onmouseover\", \"showGenomesForNode(\'\" + treeJSON[i][j][7] + \"\')\");
									newSVG.setAttributeNS(null, \"onmouseout\", \"restoreNode('\" + treeJSON[i][j][7] + \"\')\");
									newSVG.setAttributeNS(null, \"onclick\", \"keepNodeOn(\'\"+ treeJSON[i][j][7] + \"\')\");
									treeSVG.appendChild(newSVG);
								}
							}
						}

						//Adding metagenomic information to the eddges: shown as dots outside of the leaves

						//Current x-value is the number a met
						var curr_x = meta_cnt + parseInt(level)-1;
						//k = metadata variable name
						for (var k in metadata_head)
						{
							if (metadata_is_checked[k] == true)
							{
								curr_x = curr_x - 1;
								//Going through every node/branch end at the
								for (var j = 0; j < treeJSON[parseInt(level)].length; j++)
								{
									if (treeJSON[parseInt(level)][j][6] == 0)
									{
										var cx = Math.cos(parseFloat(treeJSON[parseInt(level)][j][3]) * 2 * pi) * x_dist * (levs - curr_x) + xc;
										var cy = Math.sin(parseFloat(treeJSON[parseInt(level)][j][3]) * 2 * pi) * x_dist * (levs - curr_x) + yc;
										var r = x_dist / 3;
										if (treeJSON[parseInt(level)][j][7] in node2meta)
										{
											var node_counts = node2meta[treeJSON[parseInt(level)][j][7]][k];
											var tot = treeJSON[parseInt(level)][j][4];
											var curr_p = 0;

											//Draws pie charts is utmost node < 100% of the same metadata value
											for (var i in node_counts)
											{

												var largeArc = 0; //Allows SVG to draw arc > 180 deg
												var perc = parseFloat(node_counts[i]);
												if ((perc / tot) > 0.49)
												{
													largeArc = 1;
												}

												if (perc != tot)
												{
													var x1 = Math.cos(curr_p) * r + cx;
													var y1 = Math.sin(curr_p) * r + cy;
													var x2 = Math.cos((curr_p + perc / tot) * 2 * Math.PI) * r + cx;
													var y2 = Math.sin((curr_p + perc / tot) * 2 * Math.PI) * r + cy;
													curr_p = curr_p + perc / tot;
													var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"path\");
													newSVG.setAttributeNS(null, \"id\", j + \"3\");
													newSVG.setAttributeNS(null, \"d\", \"M\" + cx+\" \" +cy+\" L \" + x1+\" \" +y1+\" A \" + r + \" \" + r + \" 0 \"+ largeArc+ \" 1 \" +x2+\" \" + y2 + \" L \" + cx + \" \" + cy);
													newSVG.setAttributeNS(null, \"stroke\", \"red\");
													newSVG.setAttributeNS(null, \"stroke-width\", \"0\");
													newSVG.setAttributeNS(null, \"fill\", metadata_head[k][i]);

													treeSVG.appendChild(newSVG);
												}
												else
												{
													var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
													newSVG.setAttributeNS(null, \"id\", j + \"3\");
													newSVG.setAttributeNS(null, \"cx\", cx);
													newSVG.setAttributeNS(null, \"cy\", cy);
													newSVG.setAttributeNS(null, \"r\", r);
													newSVG.setAttributeNS(null, \"fill\", metadata_head[k][i]);

													treeSVG.appendChild(newSVG);
												}
											}
										}
										else
										{
											var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\");
											newSVG.setAttributeNS(null, \"id\", j + \"3\");
											newSVG.setAttributeNS(null, \"cx\", cx);
											newSVG.setAttributeNS(null, \"cy\", cy);
											newSVG.setAttributeNS(null, \"r\", r);
											newSVG.setAttributeNS(null, \"fill\", \"black\");

											treeSVG.appendChild(newSVG);

										}
									}
								}
							}
						}
					}


				}
				//Make sures that the tree is drawn appropiatly and all the currently on genes
				//and nodes are kept on
				if (curr_gene != null)
				{
					add_to_tree(curr_gene);
				}
				if (nodeIsClicked != 0)
				{
					var oldNode = nodeIsClicked;
					showGenomesForNode(nodeIsClicked);
					nodeIsClicked = oldNode;
				}
			}
		}

		//
		function restoreNode(gene_id)
		{
			if (gene_id != null && gene_id in node2genome && nodeIsClicked== 0)
			{
				var insertTd = document.getElementById(\"gene_list\");
				var circ =  document.getElementsByClassName(gene_id+\"2\");
				if (circ[0] != null)
				{
					circ[0].setAttributeNS(null, \"fill\", circ[0].getAttributeNS(null, \"oldFill\"));
					var oldR = circ[0].getAttributeNS(null, \"r\");
					circ[0].setAttributeNS(null, \"r\", oldR / 2);

				}
				insertTd.innerHTML = \"\";
				var insertPd = document.getElementById(\"percDiv\");
				insertPd.innerHTML = \"\";
			}
		}

		//Shows the counts of genomes that are (a) children of a given node, (b) contain all turned on metadata variables
		// and (c) contain any selected genes/fgi in the count_genomes window on the tree page/sub page
		function showGenomesForNode(gene_id)
		{
			//gene_id is the node to be shown

			//Since only one node can be shown at a time, any node already on is turned off
			if (nodeIsClicked != 0)
			{
				var oldID = nodeIsClicked
				nodeIsClicked = 0;
				restoreNode(oldID);


			}

			//only genomes with an associated set of nodes are able to be shown
			if (gene_id != null && gene_id in node2genome)
			{
				//Getting the id of the window
				var insertTd = document.getElementById(\"gene_list\");
				//Resetting the window to blank
				insertTd.innerHTML = \"\";

				//Getting all the genomes descending from the node
				var genomeNewList = node2genome[gene_id].split(\",\");
				var genomeNewStr = \"\";
				var curr_x = 0;
				//Going through the genomes, adding one after another
				//line breaks are added at each end of a line if it exceeds 25 character
				for (var i =0; i < genomeNewList.length; i++)
				{
					if (curr_x + genomeNewList[i].length > 25)
					{
						genomeNewStr = genomeNewStr + \"<br />\";
						curr_x = 0;
					}
					genomeNewStr = genomeNewStr + genomeNewList[i] + \",\";
					curr_x = curr_x + genomeNewList[i].length + 1;
				}

				//getting the SVG circle associated with the node selected
				var circ =  document.getElementsByClassName(gene_id+\"2\");
				if (circ[0] != null)
				{
					//Changing the color to green; saving the old fill volumes to fill
					circ[0].setAttributeNS(null, \"oldFill\", circ[0].getAttributeNS(null, \"fill\"));
					circ[0].setAttributeNS(null, \"fill\", \"green\");
					//Also doubling the radius size
					var oldR = circ[0].getAttributeNS(null, \"r\");
					circ[0].setAttributeNS(null, \"r\", oldR * 2);

				}
				//Adding the genome list to the window
				insertTd.innerHTML = genomeNewStr;


				var insertPD = document.getElementById(\"percDiv\");
				insertPD.innerHTML = \"\";
				var str = \"\";

					var tot = null;
					var cnt = null;
					var circ = document.getElementsByClassName(gene_id + \"2\");
					if (circ != null)
					{
						for (var k =0 ;k < circ.length; k++)
						{
							if (circ[k].hasAttribute(\"totCnt\"))
							{
								tot = circ[k].getAttribute(\"totCnt\");
								cnt = circ[k].getAttribute(\"genomeCnt\");
							}
						}
					}
					if (tot != null)
					{
						str = str + \"Count:\" + cnt + \" (\" + parseInt(10000*cnt/tot)/100 +\")%\";
					}
					for (var i in metadata_is_checked)
					{
						if (metadata_is_checked[i])
						{
							var metaCnt = {};
							for (var j =0; j < genomeNewList.length; j++)
							{
								if (metaCnt !== undefined && metadata[genomeNewList[j]] !== undefined && metadata[genomeNewList[j]][i] in metaCnt)
								{
									metaCnt[metadata[genomeNewList[j]][i]] = metaCnt[metadata[genomeNewList[j]][i]] + 1;
								}
								else
								{
									if (metadata[genomeNewList[j]] !== undefined)
									{
										metaCnt[metadata[genomeNewList[j]][i]] = 1;
									}
								}
							}

							str = str + \"<br />\" + i + \": \";
							for (j in metaCnt)
							{
								str = str + j + \"(\" + metaCnt[j] + \"),\";
							}
							str.slice(0,-1);
						}
					}
				insertPD.innerHTML = str;

			}
		}

		//Given a list of nodes, highlights them on the tree by re-drawing the tree
		function show_nodes_on_tree(new_list)
		{
			var treeSVG = document.getElementById(\"tree_svg\");
			curr_fgi = null;

			//saves this for the future, but need to set curr_gene to null for the tree drawing
			var tmp_curr_gene = curr_gene;
			curr_gene = null;

			//redraw a blank tree
			draw_tree(tree_level, curr_type);

			//since new_list is a comma separated list of genomes (ie leaves), this makes an arrau
			var genomes = new_list.split(\",\");

			//array holding the counts of genomes
			var node_counts = [];

			//For each genome
			for (i in genomes)
			{
				if (genomes[i] in genome2node)
				{
					var node_list = genome2node[genomes[i]].split(\";\");
					for (j in node_list)
					{
						if (node_list[j] in node_counts)
						{
							node_counts[node_list[j]] = node_counts[node_list[j]] + 1;
						}
						else
						{
							node_counts[node_list[j]] =1;
						}
					}
				}
			}
			for (j in node2genome)
			{
				if (j in node_counts)
				{
					var line = document.getElementsByClassName(j + \"0\");
					if (line != null)
					{
						for (var k =0 ;k < line.length; k++)
						{
							line[k].setAttributeNS(null,\"stroke-width\", 4);
						}
					}
					var line = document.getElementsByClassName(j + \"1\");
					if (line != null)
					{
						for (var k =0 ;k < line.length; k++)
						{
							line[k].setAttributeNS(null,\"stroke-width\", 4);
						}
					}

					//get all the nodes (these are circles)
					var circ = document.getElementsByClassName(j + \"2\");
					if (circ != null)
					{
						for (var k =0 ;k < circ.length; k++)
						{
							var tot = parseInt(circ[k].getAttributeNS(null, \"cnt\"));
							if (node_counts[j] === tot)
							{
								if (j != nodeIsClicked)
								{
									circ[k].setAttributeNS(null, \"fill\", \"black\");

								}
								var old_r = parseFloat(circ[k].getAttribute(\"r\"));
								circ[k].setAttributeNS(null, \"r\", old_r * 2);
								circ[k].setAttributeNS(null, \"old_r\", old_r);
								circ[k].setAttributeNS(null, \"totCnt\", tot);
								circ[k].setAttributeNS(null, \"genomeCnt\", node_counts[j]);
							}
							else
							{
								circ[k].setAttributeNS(null, \"stroke\", circ[k].getAttribute(\"fill\"));

								circ[k].setAttributeNS(null, \"fill\", \"LightGray\");

							var cx = parseFloat(circ[k].getAttributeNS(null, \"cx\"));
								var cy = parseFloat(circ[k].getAttributeNS(null, \"cy\"));
								var r = parseFloat(circ[k].getAttributeNS(null, \"r\"))*2;
								circ[k].setAttributeNS(null, \"r\", r);
								circ[k].setAttributeNS(null, \"totCnt\", tot);
								circ[k].setAttributeNS(null, \"genomeCnt\", node_counts[j]);
								var largeArc = 0;
								if ((node_counts[j] / tot) > 0.5)
								{
									largeArc = 1;
								}
								var x1 = Math.cos(0) * r + cx;
								var y1 = Math.sin(0) * r + cy;
								var x2 = Math.cos(node_counts[j] / tot * 2 * Math.PI) * r + cx;
								var y2 = Math.sin(node_counts[j] / tot * 2 * Math.PI) * r + cy;

								var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"path\");
								newSVG.setAttributeNS(null, \"id\", j + \"3\");
								newSVG.setAttributeNS(null, \"d\", \"M\" + cx+\" \" +cy+\" L \" + x1+\" \" +y1+\" A \" + r + \" \" + r + \" 0 \"+ largeArc+ \" 1 \" +x2+\" \" + y2 + \" L \" + cx + \" \" + cy);
								newSVG.setAttributeNS(null, \"stroke\", \"red\");
								newSVG.setAttributeNS(null, \"stroke-width\", \"0\");
								newSVG.setAttributeNS(null, \"fill\", \"black\");
								newSVG.setAttributeNS(null, \"onmouseover\", circ[k].getAttributeNS(null, \"onmouseover\"));
								newSVG.setAttributeNS(null, \"onmouseout\",circ[k].getAttributeNS(null, \"onmouseout\"));
								newSVG.setAttributeNS(null, \"onclick\", circ[k].getAttributeNS(null, \"onclick\"));

								if (j == nodeIsClicked)
								{
									newSVG.setAttributeNS(null, \"fill\", \"green\");
								}
								treeSVG.appendChild(newSVG);

							}
						}
					}

				}
				else
				{
					var circ = document.getElementsByClassName(j + \"2\");
					for (var k =0 ;k < circ.length; k++)
					{

						var old_r = parseFloat(circ[k].getAttribute(\"r\"));
						circ[k].setAttributeNS(null, \"r\", old_r * 2);
						circ[k].setAttributeNS(null, \"stroke\", circ[k].getAttribute(\"fill\"));
						circ[k].setAttributeNS(null, \"fill\", \"LightGray\");
					}
				}
			}
			curr_gene = tmp_curr_gene;

		}

		//Changes the phylogenetic level shown on the tree
		//level is either - or +
		function changeTree(level)
		{
			//Move towards the leafs except if already at the bottom
			if (level == \"-\" && tree_level > 1)
			{
				//reset the level and then redraw the tree
				tree_level = tree_level -1;
				draw_tree(tree_level, curr_type);

				//Make sure that the same nodes are kept clicked/shown
				if (curr_genome != null)
				{
					add_fgi_to_tree(curr_genome);
				}
				if (nodeIsClicked != 0)
				{
					restoreNode(nodeIsClicked);
				}

				//Resetting the intra-html values of the tree level
				var p = document.getElementById(\"LevelCount\");
				p.innerHTML = tree_level;
			}

			//Move toward the root except if already there
			if (level == \"+\" && tree_level < Object.keys(treeJSON).length)
			{
				//reset the level and then redraw the tree
				tree_level = tree_level +1;
				draw_tree(tree_level, curr_type);

				//Make sure that the same nodes are kept clicked/shown
				if (curr_genome != null)
				{
					add_fgi_to_tree(curr_genome);
				}
				if (nodeIsClicked != 0)
				{
					restoreNode(nodeIsClicked);
				}

				//Resetting the intra-html values of the tree level
				var p = document.getElementById(\"LevelCount\");
				p.innerHTML = tree_level;
			}
			var genomes = document.getElementById(\"GenomesValue\").innerHTML;
			var new_list = genomes.replace(/<br>/g, \"\");
			show_nodes_on_tree(new_list);
		}

		//Upon clicking on or off the metadata values, rewrite the legend and the tree
		function change_metadata_status(evt)
		{
			var target = evt.target;
			metadata_is_checked[target.value] = target.checked;
			//redrawing legend
			change_legend();
			//redrawing trees
			draw_tree(tree_level, curr_type);

			var genomes = document.getElementById(\"GenomesValue\").innerHTML;
			var new_list = genomes.replace(/<br>/g, \"\");
			//rewriting the statistics shown in the node mouseover/click
			show_nodes_on_tree(new_list);

		}

		/*Changes the Phylogeny Metadata Legend after the Metadata sections are clicked on or off*/
		function change_legend()
		{
			//Current x-axis location in the legend
			var curr_x = 10;

			//Everytime completely remakes the legend
			var newMod = document.getElementById(\"LegendSVG\");
			newMod.innerHTML = \"\";

			//Getting the metadata values that are selected from the clicked metadata array

			for (var i in metadata_is_checked)
			{
				//Resetting x-axis variable for every Metadata variable
				curr_x = 0;
				if (metadata_is_checked[i])
				{

					var metaCheck = document.getElementById(\"select\"+i);
					var cur_y = metaCheck.offsetTop;

					//Foreach variable in the metadata, draw a circle and give the variable name in horizontal order
					//color of the circle is preset and stored in metadata_head
					//j is the string of the metadata variable name
					for (var j in metadata_head[i])
					{
						//Draw the circle showing the color and append to the legend SVG
						newModMod = document.createElementNS(\"http://www.w3.org/2000/svg\", \"circle\")
						newModMod.setAttribute(\"cx\", curr_x+5);
						newModMod.setAttribute(\"cy\", cur_y - 28);
						newModMod.setAttribute(\"r\", 4);
						newModMod.setAttribute(\"fill\", metadata_head[i][j]);
						newMod.appendChild(newModMod);

						curr_x = curr_x + 10

						newModMod = document.createElementNS(\"http://www.w3.org/2000/svg\", \"text\")
						newModMod.setAttribute(\"x\", curr_x+2);
						newModMod.setAttribute(\"y\", cur_y -22 );

						newModMod.innerHTML = j;
						newMod.appendChild(newModMod);

						//Used to calculate width of the text to advance
						var bbox = newModMod.getBBox();
						curr_x = curr_x + bbox.width;

					}


				}

			}
		}

	/*Function that decides whther to turn a tree node on or off dependening on whether a node is already clicked*/
	function keepNodeOn(gene_id)
	{
		if (nodeIsClicked == 0)
		{
			nodeIsClicked = gene_id;
		}
		else
		{
			nodeIsClicked = 0;
		}
	}

	/*Writes the Phylogeny Figure into an SVG or PNG file*/
	function saveSVGtree(fileType, id)
	{
		var event = new Event(\'build\');
		//Making a new canvas which is used to save the image to a file
		var nw = document.createElement(\"a\");

		//Loading the SVG image using the ID tree_svg
		var svg_all = document.getElementById(\"tree_svg\");
		var svgStr = \"<svg xmlns=\\\"http://www.w3.org/2000/svg\\\" xmlns:xlink=\\\"http://www.w3.org/1999/xlink\\\" height=\\\"$SVGHEIGHT\\\" width=\\\"$SVGWIDTH\\\">\"+document.getElementById(\"LegendSVG\").innerHTML+ document.getElementById(\"tree_svg\").innerHTML+\"</svg>\";
		var inSvg = new Blob([svgStr], {type: \'image/svg+xml\'});
		var svgfile = URL.createObjectURL(inSvg);

		//Getting the location to save the image to on the users computer, using the current location of the file
		var location = window.location.href;
		var curPath = location.substring(0, location.lastIndexOf(\"/\")+1);

		if (fileType==\"png\")
		{
			var newSVG = document.createElementNS(\"http://www.w3.org/2000/svg\", \"svg\");
			newSVG.setAttribute(\"xmlns:xlink\",\"http://www.w3.org/1999/xlink\");
			newSVG.setAttribute(\"height\", parseFloat(svg_all.getAttribute(\"width\")) * 2);
			newSVG.setAttribute(\"width\",2*parseFloat(svg_all.getAttribute(\"height\")) + 200);
			newSVG.innerHTML = document.getElementById(\"LegendSVG\").innerHTML+ document.getElementById(\"tree_svg\").innerHTML;
			var svgStr = new XMLSerializer().serializeToString(newSVG);

			var canvas = document.createElement(\"canvas\");
			canvas.width = parseFloat(svg_all.getAttribute(\"width\")) * $DPI;
			canvas.height = (parseFloat(svg_all.getAttribute(\"height\")) + 100)*$DPI;
			var ctx = canvas.getContext(\"2d\");
			ctx.scale($DPI, $DPI);

			var img = new Image();
			img.src = \"data:image/svg+xml;base64,\" + window.btoa(svgStr);
			
			img.onload = function()
			{
				ctx.drawImage(img, 0, 0);

				nw.setAttribute(\"download\",id + \".tree.$filenamePNG\");
				nw.setAttribute(\"href\", canvas.toDataURL(\"image/png\", 1.0));
				nw.setAttribute(\"target\", \"_blank\");

				alert(\"Saved PNG to \" + id + \".tree.$filenamePNG\");
				document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
			}
		}
		else
		{
			nw.setAttribute(\"download\", id + \".tree.$filenameSVG\");
			nw.setAttribute(\"href\", svgfile);
			nw.setAttribute(\"target\", \"_blank\");
			alert(\"Saved SVG to \" + id + \".tree.$filenameSVG\");
			document.body.appendChild(nw); nw.click(); document.body.removeChild(nw);
		}
	}
";
    close(FILE);

}

#Reading in a function file and setting the variables. Does it modularly
sub read_function_file {

    if ( !$func_file ) { die("A config file is needed. Please refer to the README"); }
    open( FILE, "<", $func_file ) or die("Cannot locate function configure file $func_file. Please check...\n");

    #Gets the different IDs
    my %map_func_ids = ( "ID" => 1, "mapfile" => 2, "ontology" => 3, "name" => 4, "color" => 5 )
      ;    #The expected variable names in the function file and where in the array to put them
    my $id  = 1;
    my $bad = 1;
    my @tmp
      ; #bad is called if (a) there isn't a start; or (b) if there isn't a recognizable variable. Prevents the information from being added into the function array
    while ( !eof(FILE) ) {

     #Only reads in the modules that start with START and ends with END and has the expected variable using the $bad variable
        my $in = <FILE>;
        chop($in);
        if ( $in =~ /END/ ) { $bad = 1; }
        if ( $in =~ /START/ ) { $bad = 0; $id++; }

        #maybe should use split....
        if ( $in =~ /([^\t]+)\t([^\t\r\n]+)/ ) {

            my $var = $1;
            my $val = $2;
            if ( !$bad ) {

                if ( !$map_func_ids{$var} ) {

                    warn
"Function configure file has unreconizable variable $var. Please refer to README. Not reading in this module...\n";
                    $id--;
                    $bad = 1;

                } else {

                    $funcInfo[ $id - 2 ]->[ $map_func_ids{$var} - 1 ] = $val;

                }

            }

        }

    }
    close(FILE);
    warn "Read in $func_file. Found ", ( $id - 1 ), " additional functional annotations...\n";

}

#Reading in a graphics config file and setting the variable
sub read_config_file {

    if ( -e $graphic_config_file ) {

        open( FILE, "<", $graphic_config_file );
        my $tmp
          ;   #getting all the variable information. changing variable name to have only uppercase and underline if neccesary
        while ( !eof(FILE) ) {

            my $v = <FILE>;
            chomp($v);
            my @n = split "\t", $v;
            $n[0] =~ s/(\s+)/\_/g;
            $tmp->{ uc( $n[0] ) } = $n[1];

        }

        #All the non-expected variables are ignored
        if ( $tmp->{BORDER_SIZE} ) { $border     = $tmp->{BORDER_SIZE}; }
        if ( $tmp->{MAX_RADIUS} )  { $max_radius = $tmp->{MAX_RADIUS}; }

        if ( $tmp->{BACKGROUND} ) { $rgb{background} = $tmp->{BACKGROUND}; }
        if ( $tmp->{fGR} )        { $rgb{fGR}        = $tmp->{fGR}; }
        if ( $tmp->{TEXT_SIZE} )  { $text_height     = $tmp->{TEXT_SIZE}; }

        if ( $tmp->{GENE_SIZE} ) { $gene_height = $tmp->{GENE_SIZE}; }

        #Using variables to calculate other specs as needed
        $layer_height = $max_radius / 12;

        $tick_height  = $gene_height / 3;
        $SVGHEIGHT    = 2 * ( $border + $max_radius );
        $SVGWIDTH     = 2 * ( $border + $max_radius );
        $thumb_size_h = $SVGHEIGHT * $thumb_dif;
        $thumb_size_w = $SVGWIDTH * $thumb_dif;
        $chng_w       = $gene_height + $thumb_size_w;
        $bw           = ( $max_radius * 2 - 2.25 * $border );
        $bodyHeight   = ( $max_radius / 4 - 30 ) . "px";

    }

}

#Quick way to get the disk image SVG, either for direct to HTML or via Javascript
sub disk_image {

    if ( $_[0] eq "default" ) {

        return "<g xmlns=\"http://www.w3.org/2000/svg\" id=\"disk_image\" ACTION>
        <path transform=\"TRANS\" style=\"fill:#434d57;fill-opacity:1;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" d=\"m 137.44242,24.150094 c 5.80716,-5.80716 20.5286,-3.244018 20.5286,-3.244018 l 619.51752,0 c 0,0 14.50156,-1.311524 19.15789,3.18591 4.05712,3.918671 3.35359,16.586114 3.35359,16.586114 l 0,664.20769 c 0,0 0.56917,10.122 -2.6273,13.33043 -3.21174,3.22375 -13.38414,2.68986 -13.38414,2.68986 l -621.66973,0 -28.98552,-29.71015 0,-640.942028 c 0,0 -2.80826,-19.18646 4.10909,-26.103808 z\" id=\"rect2826\" sodipodi:nodetypes=\"cccaccaccccs\" xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"/>
        <path style=\"fill:#949fa5;fill-opacity:1;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" d=\"m 290.34375,481.9375 0,237 350.5625,0 0,-237 -350.5625,0 z M 346,499.46875 l 85.34375,0 0,183.65625 -85.34375,0 0,-183.65625 z\" transform=\"TRANS\" id=\"rect3607\"/>
        <rect style=\"fill:#eeeeee;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" id=\"rect3601\" width=\"523.46179\" height=\"392.23584\" x=\"209.09631\" y=\"34.609047\" transform=\"TRANS\" rx=\"20\" ry=\"20\"/>
        <rect style=\"fill:#333333;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" id=\"rect3603\" width=\"37.342163\" height=\"28.197144\" x=\"157.75159\" y=\"70.411377\" transform=\"TRANS\"/>
        <rect transform=\"TRANS\" y=\"378.67908\" x=\"749.12952\" height=\"28.197144\" width=\"37.342163\" id=\"rect3605\" style=\"fill:#333333;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"/>
		<text transform=\"TRANS\" x=\"450\" y=\"200\" text-anchor=\"middle\" dominant-baseline=\"middle\" stroke-width=\"20\" stroke=\"#0000dd\" font-size=\"FS\">TEXT</text></g>";

    }
    if ( $_[0] eq "js" ) {

        return
"<g xmlns=\\\"http://www.w3.org/2000/svg\\\" id=\\\"disk_image\\\" ACTION><path transform=\\\"TRANS\\\" style=\\\"fill:#434d57;fill-opacity:1;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\\\" d=\\\"m 137.44242,24.150094 c 5.80716,-5.80716 20.5286,-3.244018 20.5286,-3.244018 l 619.51752,0 c 0,0 14.50156,-1.311524 19.15789,3.18591 4.05712,3.918671 3.35359,16.586114 3.35359,16.586114 l 0,664.20769 c 0,0 0.56917,10.122 -2.6273,13.33043 -3.21174,3.22375 -13.38414,2.68986 -13.38414,2.68986 l -621.66973,0 -28.98552,-29.71015 0,-640.942028 c 0,0 -2.80826,-19.18646 4.10909,-26.103808 z\\\" id=\\\"rect2826\\\" sodipodi:nodetypes=\\\"cccaccaccccs\\\" xmlns:sodipodi=\\\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\\\"/><path style=\\\"fill:#949fa5;fill-opacity:1;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\\\" d=\\\"m 290.34375,481.9375 0,237 350.5625,0 0,-237 -350.5625,0 z M 346,499.46875 l 85.34375,0 0,183.65625 -85.34375,0 0,-183.65625 z\\\" transform=\\\"TRANS\\\" id=\\\"rect3607\\\"/><rect style=\\\"fill:#eeeeee;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\\\" id=\\\"rect3601\\\" width=\\\"523.46179\\\" height=\\\"392.23584\\\" x=\\\"209.09631\\\" y=\\\"34.609047\\\" transform=\\\"TRANS\\\" rx=\\\"20\\\" ry=\\\"20\\\"/><rect style=\\\"fill:#333333;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\\\" id=\\\"rect3603\\\" width=\\\"37.342163\\\" height=\\\"28.197144\\\" x=\\\"157.75159\\\" y=\\\"70.411377\\\" transform=\\\"TRANS\\\"/><rect transform=\\\"TRANS\\\" y=\\\"378.67908\\\" x=\\\"749.12952\\\" height=\\\"28.197144\\\" width=\\\"37.342163\\\" id=\\\"rect3605\\\" style=\\\"fill:#333333;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\\\"/><text transform=\\\"TRANS\\\" x=\\\"450\\\" y=\\\"200\\\" text-anchor=\\\"middle\\\" dominant-baseline=\\\"middle\\\" stroke-width=\\\"20\\\" stroke=\\\"#0000dd\\\" font-size=\\\"FS\\\">TEXT</text></g>";

    }

}

#Running a color wheel function to get a certain set of opposing color values reported as rgb
sub get_rgb_string {

    my $in = $_[0];
    my @arr = ( 0, 0, 0 );
    if ( $in >= 0 && $in < 60 ) {

        $arr[0] = 255;
        $arr[1] = $in * 255 / 60;

    }
    if ( $in >= 60 && $in < 120 ) {

        $arr[1] = 255;
        $arr[0] = ( 120 - $in ) * 255 / 60;

    }

    if ( $in >= 120 && $in < 180 ) {

        $arr[1] = 255;
        $arr[2] = ( $in - 120 ) * 255 / 60;

    }

    if ( $in >= 180 && $in < 240 ) {

        $arr[2] = 255;
        $arr[1] = ( 240 - $in ) * 255 / 60;

    }
    if ( $in >= 240 && $in < 300 ) {

        $arr[2] = 255;
        $arr[0] = ( $in - 240 ) * 255 / 60;

    }

    if ( $in >= 300 && $in < 360 ) {

        $arr[0] = 255;
        $arr[2] = ( 360 - $in ) * 255 / 60;

    }
    if ( $in == 360 ) { $arr[0] = 255; }
    return ( sprintf( "#%.2x%.2x%.2x", $arr[0], $arr[1], $arr[2] ) );

}

#Loading the tree from bioperl into the bioperl format
sub read_tree {

    my $outfile = $_[0];
    my $treeio = Bio::TreeIO->new( -format => 'newick', -file => $outfile );
    if ( !$treeio ) { print "Cannot read in tree from $outfile.... removing tree but continuing\n\n"; return (0); }
    my $tree = $treeio->next_tree;
    if ( !$tree ) { print "Cannot read in tree $tree from $outfile... removing tree but continuing\n\n"; return (0); }
    return $tree;

}

#Making an array from a BioPerl tree format. Is recursive
#Makes a hash with the keys the BioPerl node ID. Secondary hash keys are:
#		depth (the level of the node: 1 = leaf)
#		descendent_list (string containing list of leaves descending from this node)
#		total descendents (the number of leaves descending from this node)
#		reverse depth (the level of the node where 1 = root)
#		parent (hash key of the parent node)
#		child (count of the the number of children)
#		child_ids (semi-colon seperated list of hash keys for the direct descendents of the node)

sub tree_to_array {

    my $tree = $_[0];                 #pointer to root ID
    my $arr  = $_[1];                 #pointer to hash
    my $id   = $tree->internal_id;    #Getting the ID

    #Adding the name to the tree (this is primarily leave ids)
    if ( $tree->id ) {

        $arr->{$id}->{name} = $tree->id;

    }

    #Setting the leaf variables: Depth and des
    if ( $tree->is_Leaf ) {

        $arr->{$id}->{depth}             = 1;
        $arr->{$id}->{total_descendents} = 1;
        $arr->{$id}->{descendent_list}   = "<" . $tree->id . ">";

    }

    #Check for parent. If there isn't a parent, then the node is a root
    if ( $arr->{$id}->{parent} ) {

        $arr->{$id}->{reverse_depth} = $arr->{ $arr->{$id}->{parent} }->{reverse_depth} + 1;

    } else {

        $arr->{$id}->{reverse_depth} = 1;

    }

    #Looking at each child of the nodes
    foreach my $node ( $tree->each_Descendent ) {

        my $id2 = $node->internal_id;
        if ( !$arr->{$id2} ) {

            $arr->{$id2}->{parent} = $id;
            $arr->{$id}->{child}++;
            $arr->{$id}->{child_ids} .= "$id2;";
            $arr = tree_to_array( $node, $arr );

        }
        if ( $arr->{$id}->{depth} ) {

            $arr->{$id}->{depth} = max( $arr->{$id}->{depth}, $arr->{$id2}->{depth} + 1 );

        } else {

            $arr->{$id}->{depth} = $arr->{$id2}->{depth} + 1;

        }
        $arr->{$id}->{total_descendents} += $arr->{$id2}->{total_descendents};
        $arr->{$id}->{descendent_list} .= $arr->{$id2}->{descendent_list};

    }

    return ($arr);

}

#Makes an array out of the phylogeny to be written into a json file
#out is an array of arrays with the first key is the level and the second are the different parts (nodes & branches)
#row is [x1, y1, x2, y2, list_of_children_genome, list_of_all_genomes, type (0 = vertical branch, 1 = horizontal branch, 2 = node), node id]
sub make_tree_array {

    my $arr = $_[0];
    my $loc;

    #out is an array of arrays with the first key is the level and the second are the different parts (nodes & branches)
    my $out;

    #sorting the tree array by depth
    my @s = sort { $arr->{$b}->{depth} <=> $arr->{$a}->{depth} } keys(%$arr);

    #Getting the array key of the root
    my $a = $s[0];
    $loc->{$a}->{mn_y}  = 0;
    $loc->{$a}->{mx_y}  = 1;
    $loc->{$a}->{st_x}  = $arr->{$a}->{depth} + 1;
    $loc->{$a}->{end_x} = $arr->{$a}->{depth};

    #The total number of genomes/leaves of the tree
    my $tot = $arr->{$a}->{total_descendents};

    #@list is array of the nodes. Starting by pushing the root onto the list
    my @list;
    push @list, $a;

    #max_x = depth
    my $max_x = $arr->{$a}->{depth} + 1;

    #putting in the root node
    $out->[ $max_x - 1 ]->[0] = [ $max_x, 0.5, $max_x - 1, 0.5, $tot, $tot, 0, $a ];

    #Going through the list (it is iteratively added to), and adding the nodes/branches/etc into the out array
    while ( scalar(@list) > 0 ) {

        my $a     = $list[0];
        my $cur_y = $arr->{$a}->{depth};    #starting depth
        my @n;                              #array with children nodes
        if ( $arr->{$a}->{child_ids} ) { @n = split ";", $arr->{$a}->{child_ids}; }
        my $mn_y = $loc->{$a}->{mn_y};                                 #starting horizontal point
        my $mn_x = $loc->{$a}->{end_x};                                #ending of the horizontal point
        my $y1   = ( $loc->{$a}->{mn_y} + $loc->{$a}->{mx_y} ) / 2;    #y1 is the middle horizontal point

        #Goes through all the vertical levels to add the branch one at a time
        for ( my $i = $loc->{$a}->{end_x} ; $i < $loc->{$a}->{st_x} ; $i++ ) {

            if ( $out->[$i] ) {

                push @{ $out->[$i] }, [ $i, $y1, ( $i + 1 ), $y1, $arr->{$a}->{total_descendents}, $tot, 0, $a ];

            } else {

                #First row of the level
                $out->[$i]->[0] = [ $i, $y1, ( $i + 1 ), $y1, $arr->{$a}->{total_descendents}, $tot, 0, $a ];

            }

        }

        #adding node
        if ( $out->[ $loc->{$a}->{end_x} ] ) {

            push @{ $out->[ $loc->{$a}->{end_x} ] },
              [ $loc->{$a}->{end_x}, $y1, $loc->{$a}->{end_x}, $y1, $arr->{$a}->{total_descendents}, $tot, 2, $a ];

        } else {

            $out->[ $loc->{$a}->{end_x} ]->[0] =
              [ $loc->{$a}->{end_x}, $y1, $loc->{$a}->{end_x}, $y1, $arr->{$a}->{total_descendents}, $tot, 2, $a ];

        }
        my $md_y;    #
        my $st_y = $y1;    #

        #Going through all the children of the node
        foreach my $n (@n) {

            push @list, $n;    #Pushing the children onto list

            $loc->{$n}->{end_x} = $arr->{$n}->{depth};    #vartical beginning
            $loc->{$n}->{st_x}  = $mn_x;                  # vertical endnig
            my $tmp_y = $arr->{$n}->{total_descendents} / $tot;    #Percentage of genomes...
            $loc->{$n}->{mn_y} = $mn_y;                            #horizontal start

            my $i = $mn_x;                                         #level
            if ( $i > 1 ) {

                #adding horizonal level
                if ( $out->[ $i - 1 ] ) {

                    push @{ $out->[ $i - 1 ] }, [ $i, $y1, $i, $mn_y + $tmp_y / 2, -1, $tot, 1, $n ];

                } else {

                    #initial
                    $out->[ $i - 1 ]->[0] = [ $i, $y1, $i, $mn_y + $tmp_y / 2, -1, $tot, 1, $n ];

                }

            }

            $st_y = $mn_y + $tmp_y / 2;
            $mn_y += $tmp_y;
            $loc->{$n}->{mx_y} = $mn_y;

        }

        #popping the list
        shift @list;

    }
    return ($out);

}

#Translates the Bio::Tree style phylogeny into a json file that can be read by the javascript
sub make_tree_json {

    #Loading up the Bio module (treeIO)
    if ( eval { require Bio::TreeIO; 1; } ) {

        Bio::TreeIO->import();

    } else {

        warn "Cannot find BioPerl TreeIO module to load the phylogeny. Ignoring phylogeny..\n\n";
        return (0);

    }

    #Reading in the tree
    my $tree = read_tree( $_[0] );
    if ( $tree == 0 ) { return (0); }

    my $root = $tree->get_root_node;

    #Need to pre-set array then loads tree into the array
    my $arr;
    $arr = tree_to_array( $root, $arr );

    #Makes an
    my $tree_array = make_tree_array($arr);

    #counting the number of tree levels
    $tree_level = scalar(@$tree_array);

    #list out is the string giving the json hash with the genome as key and ancestral node as the value
    my $list_out = "var gene2node = \'{";

    #gets the list of leaf nodes in an array. This is a bioperl function
    my @leafs = $tree->get_leaf_nodes;

    my %leaf_2_node;    #matches leaves to the ancestral nodes
    my %leafs;          #id
    foreach my $a (@leafs) { $leafs{ $a->id } = 1; }    #matches node id as a leaf in hash
    foreach my $leaf (@leafs) {

        my $gen = $leaf->id;                            #genome id
        my $tst = "<$gen>";                             # list of genomes
        $list_out .= "\"$gen\":\"";                     #name of
        foreach my $node ( keys(%$arr) ) {

            if ( $arr->{$node}->{descendent_list} =~ /$tst/ ) {

                $list_out .= "$node;";                  #
                $leaf_2_node{$gen} .= "$node;";         #leaf to the different nodes

            }

        }
        $list_out .= "\",";

    }
    chop($list_out);
    $list_out .= "}\';";                                #json is finished

    #if available, add in metadata
    #%meta_data, $meta_out is the json string for the metadata
    #
    my ( %meta_data, $meta_out, %meta_vars, $meta_cnts, %meta_node_cnts );
    if ($meta_data_file) {

        while ( $meta_data_file =~ /([^\,]+)/g ) {

            my $meta_file_in = $1;
            if ( -e $meta_file_in ) {

                warn $meta_file_in, "\n";
                open( FILE, "<", $meta_file_in );
                while ( !eof(FILE) ) {

                    my $v = <FILE>;
                    chomp($v);
                    my @n;
                    my $c;
                    while ( $v =~ /(\S+)/g ) { $n[ $c++ ] = $1; }
                    if ( scalar(@n) == 3 && $leafs{ $n[0] } ) {

                        $meta_vars{ $n[1] }->{ $n[2] } = 1;
                        $meta_data{ $n[0] }->{ $n[1] } = $n[2];
                        my @nodes = split ";", $leaf_2_node{ $n[0] };
                        foreach my $node (@nodes) {

                            $meta_node_cnts{$node}->{ $n[1] }->{ $n[2] }++;

                        }

                    }

                }

            }

        }
        $meta_out .= "var meta_data = \'{";
        foreach my $a ( keys(%meta_data) ) {

            $meta_out .= "\"$a\":{";
            foreach my $b ( keys( %{ $meta_data{$a} } ) ) {

                $meta_out .= "\"$b\":\"" . $meta_data{$a}->{$b} . "\",";

            }
            chop($meta_out);
            $meta_out .= "},";

        }
        chop($meta_out);
        $meta_out .= "}\';\n";
        $meta_out .= "var meta_node_data = \'{";
        foreach my $a ( keys(%meta_node_cnts) ) {

            $meta_out .= "\"$a\":{";
            foreach my $b ( keys( %{ $meta_node_cnts{$a} } ) ) {

                $meta_out .= "\"$b\":{";
                foreach my $c ( keys( %{ $meta_node_cnts{$a}->{$b} } ) ) {

                    $meta_out .= "\"$c\":\"" . $meta_node_cnts{$a}->{$b}->{$c} . "\",";

                }
                chop($meta_out);
                $meta_out .= "},";

            }
            chop($meta_out);
            $meta_out .= "},";

        }
        chop($meta_out);
        $meta_out .= "}\';\n";
		$meta_cnts = 0;
        foreach my $a ( keys(%meta_vars) ) {

            foreach my $b ( keys( %{ $meta_vars{$a} } ) ) {

                $meta_cnts++;

            }

        }
		my $meta_col;
		if ($meta_cnts > 0)
		{
			$meta_col = 360 / ($meta_cnts);    #getting the colors for the metadata variables using get_rgb_string
        }
		$meta_cnts = 0;
        $meta_out .= "var meta_head_var = \'{";
        foreach my $a ( keys(%meta_vars) ) {

            $meta_out .= "\"$a\":{";
            foreach my $b ( keys( %{ $meta_vars{$a} } ) ) {

                $meta_out .= "\"$b\":\"" . get_rgb_string( $meta_col * $meta_cnts ) . "\",";
                $meta_types{$a}->{$b} = $meta_col * $meta_cnts;
                $meta_cnts++;

            }
            chop($meta_out);
            $meta_out .= "},";

        }
        chop($meta_out);
        $meta_out .= "}\';\n";

    }

    #Writing the tree and metadata information to the tree json file
    open( TREE_JSON, ">", "$out_dir/$output_id/json/tree.json" );

    #tree_format is the tree array in JSON format for the javascript to read and draw
    my $tree_out = "var tree_format = \'{";
    for ( my $i = 0 ; $i < scalar(@$tree_array) ; $i++ ) {

        if ( $tree_array->[$i] ) {

            $tree_out .= "\"$i\":[";

            for ( my $j = 0 ; $j < scalar( @{ $tree_array->[$i] } ) ; $j++ ) {

                $tree_out .= "[";
                for ( my $k = 0 ; $k < scalar( @{ $tree_array->[$i]->[$j] } ) ; $k++ ) {

                    $tree_out .= $tree_array->[$i]->[$j]->[$k] . ",";

                }
                chop($tree_out);
                $tree_out .= "],";

            }

            chop($tree_out);
            $tree_out .= "],";

        }

    }
    chop($tree_out);
    $tree_out .= "}\';";

    #Printing the tree array and the other tree json strings
    print TREE_JSON $tree_out, "\n", $list_out, "\n";

    #printing out the metadata json information when available
    if ($meta_out) { print TREE_JSON $meta_out, "\n"; }
    close(TREE_JSON);

}

#HTML script for the loading and waiting screens
sub make_loading_screen() {

    my $out;
    $out = sprintf(
"<style>.shifting {\n content:attr(data-content); position: absolute; overflow: hidden; color: red; max-width: 7em; animation: loading 3s linear infinite;}\n .loader {\n border: 16px solid #f3f3f3; border-top: 16px solid #3498db;  border-radius: 50%s; width: %fpx; height: %fpx;  animation: spin 2s linear infinite; }\n \@keyframes spin { 0%s { transform: rotate(0deg); } 100%s { transform: rotate(360deg); }}\n \@keyframes loading { 0%s {max-width:0; }}\n </style>",
        "%", 0, 0, "%", "%", "%" );
    $out .= sprintf(
"<style>.block1{\n animation: block1 4s linear infinite;}\n .block2{\n animation: block2 4s linear infinite;}\n .block3{\n animation: block3 4s linear infinite;}\n .block4{\n animation: block4 4s linear infinite;}\n"
    );
    $out .= sprintf(".twirl {\n animation: flip 16s linear infinite;}\n");
    $out .= sprintf(
"\@keyframes block1{ 0%s {max-width:%fpx; left:%fpx;} 12.5%s, 50%s {max-width:%fpx; left:%fpx;} 62.5%s {max-width:%fpx; left:%fpx;} 100%s {max-width:%fpx; left:%fpx;}} \n",
        "%", 0, ( $border + $max_radius - $border - 5 ),
        "%", "%", $border, ( $border + $max_radius - $border - 5 ),
        "%", 0, ( $border + $max_radius - 5 ),
        "%", 0, ( $border + $max_radius - $border - 5 )
    );
    $out .= sprintf(
"\@keyframes block2{ 0%s {max-height:%fpx; top:%fpx;} 12.5%s {max-height:%fpx; top:%fpx;} 25%s, 62.5%s {max-height:%fpx; top:%fpx;} 75%s {max-height:%fpx; top:%fpx;} 100%s {max-height:%fpx; top:%fpx;}} \n",
        "%",                                     0,
        ( $border + $max_radius - $border - 5 ), "%",
        0, ( $border + $max_radius - $border - 5 ),
        "%", "%",
        $border, ( $border + $max_radius - $border - 5 ),
        "%",                           0,
        ( $border + $max_radius - 5 ), "%",
        0, ( $border + $max_radius - $border - 5 )
    );
    $out .= sprintf(
"\@keyframes block3{ 0%s {max-width:%fpx; left:%fpx;} 25%s {max-width:%fpx; left:%fpx;} 37.5%s, 75%s {max-width:%fpx; left:%fpx;}  87.5%s {max-width:%fpx; left:%fpx;} 100%s {max-width:%fpx; left:%fpx;}} \n",
        "%",                                     0,
        ( $border + $max_radius + $border + 5 ), "%",
        0, ( $border + $max_radius + $border + 5 ),
        "%", "%",
        $border, ( $border + $max_radius + 5 ),
        "%",                           0,
        ( $border + $max_radius + 5 ), "%",
        0, ( $border + $max_radius + $border + 5 )
    );
    $out .= sprintf(
"\@keyframes block4{ 0%s {max-height:%fpx; top:%fpx;} 37.5%s {max-height:%fpx; top:%fpx;} 50%s, 87.5%s {max-height:%fpx; top:%fpx;} 100%s {max-height:%fpx; top:%fpx;}} \n",
        "%", 0, ( $border + $max_radius + $border + 5 ),
        "%", 0, ( $border + $max_radius + $border + 5 ),
        "%", "%", $border, ( $border + $max_radius + 5 ),
        "%", 0, ( $border + $max_radius + 5 )
    );
    $out .= sprintf( "\@keyframes flip {100%s {transform:rotate(360deg); }}", "%" );
    $out .= sprintf("</style>\n");
    $out .= sprintf(
"<div id='Start' style=\"visibility:visible; top:0px; left:0px; height:%fpx; width:%fpx; position: absolute;\"></div>",
        ( $border + $max_radius ) * 2,
        ( $border + $max_radius ) * 2
    );

    $out .= sprintf(
"<div id='Loading' class=\'twirl\' style=\"visibility:visible; top:%fpx; left:%fpx; height:%fpx; width:%fpx; position: fixed;\">",
        0, 0,
        2 * ( $border + $max_radius ),
        2 * ( $border + $max_radius )
    );
    $out .= sprintf(
"<div id='block1' class='block1' style='position:absolute; background:#F99D31; top:%fpx; left:%fpx; height:%fpx; width:%fpx;'></div>",
        ( $border + $max_radius - $border - 5 ),
        ( $border + $max_radius - $border - 5 ),
        $border, $border
    );
    $out .= sprintf(
"<div id='block2' class='block2' style='position:absolute; background:#6DB33F; top:%fpx; left:%fpx; height:%fpx; width:%fpx;'></div>",
        ( $border + $max_radius - $border - 5 ),
        ( $border + $max_radius + 5 ),
        $border, $border
    );
    $out .= sprintf(
"<div id='block3' class='block3' style='position:absolute; background:#E31B23; top:%fpx; left:%fpx; height:%fpx; width:%fpx;'></div>",
        ( $border + $max_radius + 5 ),
        ( $border + $max_radius + 5 ),
        $border, $border
    );
    $out .= sprintf(
"<div id='block4' class='block4' style='position:absolute; background:#00A4E4; top:%fpx; left:%fpx; height:%fpx; width:%fpx;'></div></div>",
        ( $border + $max_radius + 5 ),
        ( $border + $max_radius - $border - 5 ),
        $border, $border
    );
    return ($out);

}

#Makes the directory where the program writes the html, javascript and json files
sub make_output_dirs {

    my $curr_dir = $out_dir;
    my @add_to_dir;            #array used to store sub-directories that are needed to be made
    while ( !-e $curr_dir )    #finding all the sub-directories that need to be made
    {

        $curr_dir =~ /\A(.+)\/([^\/]+)([\/]*)\Z/;
        push @add_to_dir, $2;
        $curr_dir = $1 . "\/";

    }
    while ( scalar(@add_to_dir) > 0 ) {

        my $dir_x = pop @add_to_dir;
        mkdir( $curr_dir . "/" . $dir_x );
        $curr_dir .= $dir_x;

    }

    #The following sets up the output directory. The structure is shown in the manual
    if ( !-e "$out_dir/$output_id/" ) {

        mkdir("$out_dir/$output_id") or die("Cannot make directory....\n\n");
        warn "Making directory $out_dir/$output_id\n\n";

    }
    if ( !-e "$out_dir/$output_id/fgi/" ) {

        mkdir("$out_dir/$output_id/fgi/") or die("Cannot make directory....\n\n");
        warn "Making directory $out_dir/$output_id/fgi\n\n";

    }
    if ( !-e "$out_dir/$output_id/core/" ) {

        mkdir("$out_dir/$output_id/core/") or die("Cannot make directory....\n\n");
        warn "Making directory $out_dir/$output_id/core\n\n";

    }
    if ( !-e "$out_dir/$output_id/json/" ) {

        mkdir("$out_dir/$output_id/json/") or die("Cannot make directory....\n\n");
        warn "Making directory $out_dir/$output_id/core\n\n";

    }
    if ( !-e "$out_dir/$output_id/scripts/" ) {

        mkdir("$out_dir/$output_id/scripts/") or die("Cannot make directory....\n\n");
        warn "Making directory $out_dir/$output_id/core\n\n";

    }

}

#################################################################################################################
#
# Sub Reutines
#
# make_main_figure: makes the central pan-chromosome image.
#	Passed: None
#	Returns: None
#	Called by: init()
#	makes the SVG images for the main page by running through
#
# make_fgi_page_svg:
#	Passed: 	$file => string with the name to the file
#				$data => pointer to hash with the gene locations, order, etc
#				$ont => pointer to hash with the functional information
#				$out => pointer to hash with "leveled" svg that makes the images (starts with empty)
#	Returns: 	$out => pointer to hash with SVG images added (written to $file file)
#	Called By:	make_main_figure()
#	Makes the full SVG image for the fGI to be written to a new file in make_main_figure(). Probably could just
#
# make_fgi_region_svg
#	Passed: 	$file => string with the file name
#				$data => pointer to hash with the gene locations, order, size, etc
#				$ont => pointer to hash with the functional information
#				$out => pointer to hash with "leveled" svg that makes the images
#	Returns: 	$out => same from passed with SVG images added
#	Called By:	make_main_figure()
#	Makes the preview SVG image for the fGI to be added to
#
# make_core_region_svg
#	Passed: 	$file => string with the file,
#				$data => ppointer to hash with the gene location
#				$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: 	$out => pointer to hash with SVG images added
#	Called By:	make_main_figure()
#	Makes the full SVG image for the Core to be written in the page
#
# make_core_region_page
#	Passed: 	$file => string with the file,
#				$data => ppointer to hash with the gene location
#				$ont => pointer to hash with the functional information
#				$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: 	$out => pointer to hash with SVG images added
#	Called By:	make_main_figure()
#	Makes the full SVG image for the fGI to be written in the page
#
# svg_draw_plasmid_seg
#	Passed:		$depth => the svg "level": lower values are drawn later
#				$d1 => integer with the starting position of the gene
#				$d3 => integer with the ending position of the gene
#				$l1 => bottom height
#				$l2 => top height
#				$c => string containing the color of the plasmid segment
#				$href => string contining Plasmid ID
#				$name => Gene name
#				$id => string containing the cluster ID
#				$type => string with Gene, fGI or core cluster. Should be gene
#				$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns:	$out => with plasmid added to SVG in key level
#	Called By: make_main_figure()
#	Adds gene segment for a circular fgi (plasmid) to the SVG main page

# svg_tick
#	Passed: $x => x-axis location
#			$l => layer value
#			$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: $out => with tick added through svg_spoke()
#	Called By: make_main_figure(), place_svg_ticks(), place_svg_ticks_arc(), place_svg_ticks_multi()
#	Draws a single tick on the outside of the chromosomes

# svg_draw_arc_seg
#	Passed:	$depth => the svg "level": lower values are drawn later
#			$d1 => integer with the starting position of the gene
#			$d3 => integer with the ending position of the gene
#			$l1 => double with bottom height
#			$l2 => double with top height
#			$c => string with rgb color
#			$href => string containing location of the svg to load on click
#			$name => string with long name shown on mouse over
#			$id => string with summary name
#			$type => string showing Region or Gene
#			$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: $out => with the segment drawn in
#	Called By: make_main_figure()
#	Draws a segment on the chromosome from $d1 to $d3. It is able to account for arcs > 180 degrees

# svg_draw_arc_full
#	Passed:	$depth => the svg "level": lower values are drawn later
#			$l1 => double with bottom height
#			$l2 => double with top height
#			$c => string with rgb color
#			$href => string with the link to go to on click
#			$name => string with long name shown on mouse over
#			$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: $out => with the segment drawn in
#	Called By: make_main_figure()
#	Draws a full circle

# place_svg_ticks
#	Passed:	$increment => number with the distance between the ticks
#			$size => the size of the chrosomoe to draw it in
#			$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: $out => with the ticks drawn in from zero to size
#
#	Draws ticks at all the points from 0 to the size where size is the full length of the arc. Draws over 360 degrees

# place_svg_ticks_arc
#	Passed:	$increment => number with the distance between the ticks
#			$size => the size of the arc to draw the ticks
#			$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: $out => with the ticks drawn in from zero to size
#
#	Draws ticks at all the points from 0 to the size. Needs $seq_len as the full length

# place_svg_ticks_multi
#	Passed:	$increment => number with the distance between the ticks
#			$multi => pointer to hash with chromosome lengths
#			$out => pointer to hash with SVG that makes the images with the keys as the levels
#	Returns: $out => with the ticks drawn in from zero to size for all of the
#
#	Draws ticks at all the points from 0 to the size for a group of chromosomes all using increment to out

# svg_spoke
#	Passed: $depth => level of svg to add the stroke: not used
#			$d: x level used in the screen trans to get x and y levels
#			$l1: y start to get x and y in screen_trans
#			$l2: y end to get x and y in screen_trans
#			$c: color (not used)
#			$width: width of line (not used)
#			$href: id (not used)
#			$name: name (not used)
#			$out: svg string start that the svg path is added to
#	Returns: $out => return svg with the stroke path added
#
#	Draws a svg path of a spoke (line on an angle) to $out. Used for drawing tick

# convert_layer
#	Passed: $l
#	Returns: the height
#
#	Converts to

# screen_x_trans
#	Passed: $p => location along circle
#			$r => distance from center
#	Returns: returns x position on the circle in the total screen
#
# uses an angle ()

# screen_y_trans
#	Passed: $p => location along circle
#			$l => distance from center
#	Returns: turns y position on the circle in the total screen
#

# screen_x_trans_ratio
#	Passed: $p => location along circle
#			$r => distance from center
#			$sz => size of the circle
#	Returns: returns x position on the circle using the size passed in
#

# screen_y_trans_ratio
#	Passed: $p => location along circle
#			$r => distance from center
#			$sz => size of the circle
#	Returns: returns y position on the circle using the size passed in

# rad2deg
#	Passed: $r => radion
#	Returns: float with the value in degree
#

# deg2rad
#	Passed: $d => degree
#	Returns: float with the value in radons
#

# get_center_x
#	Passed: NONE
#	Returns: the x-positon of the center of the figure
#

# get_center_y
#	Passed: NONE
#	Returns: the y-positon of the center of the figure
#

# commas
#	Passed: number
#	Returns: string with commas in the appropiate place
#

# svg_init_image
#	Passed: None
#	Returns: string with svg header for html
#

# max
#	Passed: list of variables
#	Returns: the highest value of this
#

# min
#	Passed: list of variables
#	Returns: the lowest value of this
#

# split_id
#	Passed: $word => word used to split into multiple rows as needed
#			$x => x location
#			$len => max length allowed before making a new row
#	Returns: string with tspans with the words insides
#

# make_main_javascript
#	Passed: none
#	Returns: none
#	Writes the javascript file used by the main page in the output directory

# load_functions
#	Passed: none
#	Returns: none
#	Reads the function files passed to the program via the function config file. Writes information to global variables

# get_gene_info
#	Passed: none
#	Returns: none
#	Reads in the information from fasta files and the PanACEA Flatfile
#	Writes information to global variables

# fgi_javascript
#	Passed: none
#	Returns: none
#	Writes the javascript file used by the Core and FGI pages in the output directory

# write_geneFile_javascript
#	Passed: none
#	Returns: none
#	Writes the javascript file used by the individual gene pages in the output directory

# read_function_file
#	Passed: none
#	Returns: none
#	Reads the function config file passed to the program. Writes information to global variables

# read_config_file
#	Passed: none
#	Returns: none
#
#	Reads the function config file passed to the program. Writes information to global variables

# disk_image
#	Passed: $string => either passed "default" for an html or "js" for a javascript
#	Returns: string with disk image written in svg
#

# get_rgb_string
#	Passed: a number from 0 to 360
#	Returns: an rgb string showing the color on the colorwheel with 0/360 = red
#

# make_tree_json
#	Passed: $_[0] => string with file name
#	Returns: none
#
#	Reads in tree and metadata files and writes json file containing all the information that is
#	readable by the html files

# read_tree
#	Passed: $outfile => string with filename
#	Returns: BioPerl tree hash
#

# tree_to_array
#	Passed: $tree => pointer to root
#			$arr => pointer to array to return
#	Returns: $arr => hash with the keys as nodes
#

# make_tree_array
#	Passed: $_[0] => hash with the keys as nodes as tree
#	Returns: array of arrays with first key as the tree level and the second key
#	This is pass along to tree json and then used to write the json file

# make_loading screen
#	Passed: None
#	Returns: string with HTML script for the loading screen
#

# make_output_dirs
#	Passed: None ($out_dir is passed along)
#	Returns: none
#
#	Makes the directories to write the files
#
##########################################################################################
