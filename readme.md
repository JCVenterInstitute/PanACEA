# PanACEA (Pan-genome Atlas and Chromosome Explorer and Analyzer)
PanACEA is a suite of PERL scripts that allows users to create an interconnected set of 
html, javascript and json files that allows for visualization of prokaryotic pan-chromosomes,
including core and variable regions. Users can also include additional information such as
gene annotation, metadata annotation, and genome phylogenies. Detailed information about the PanACEA
suite and usages can be found in the PDF manual included on the site.
 
Also included in the Github is a set of exemplar data from the PanOCT run on Enterobactor hormaechei
from *Chavda et al 2016* and a shell script which can generate the PanACEA 
output files from the examples, and the gzipped output from a PanACEA run using the same data. 

## Included PERL scripts

### make_panacea.pl
make_panacea.pl generates an interconnected set of html, javascript and json files
from a PanACEA flat file generated from a possible suite of programs, including PanOCT,
that allows users to explore multiple related prokaryotic genomes through their
pan-chromosomes, including core and flexible regions. While the basic version only shows
the location and relative position of the core and flexible regions, the user can also
add information on the gene cluster annotation and alignment as well as the genomic
annotation and phylogeny to obtain a more complete view of the relationship between the
genomes, their function and any regions with their associated genes that might drive
this relationship. The end result is a central html file that can be opened in many web
browsers.

**running make_panacea.pl**

	make_panacea.pl [options]
	Options:
		-i	--input 	Input PanACEA Flat File. Required.
		-o	--output 	Output Directory for all the files. Default is current directory.
		-n	--name		Header name of HTML and SVG files. Default is \"PanHTML\"
		-d	--dir_name	Root directory of the HTML for use in the browser. Default is current directory.
		-f	--function	File containing information about additional gene cluster annotation.
		-a	--fasta		Either directory containing the cluster multiple alignment files or a single fasta file. Multiple inputs could be added in a comma separated list
		-t	--tree		File containing the phylogeny of the genomes in newick format. (Requires BioPerl to use)
		-m	--metafiles	File containing group metadata for the genomes. Only used in connection with tree.
		-g	--graphic	File containing the graphics configure file to change the size of the whole image and parts of the image

### make_panacea_flatfile.pl
make_panacea_flatfile.pl generates the panchromosome PanACEA flat file from pangenome output 
from  various programs, including PanOCT. This file is required to run the PanACEA 
visualizer.

**running make_panacea.pl**

	make_panacea_flatfile.pl [options]
	Options:
		-i	--input 	Input directory of the Pangenome information. Required.
		-o	--output 	Output file. Default is PanACEA.[date].txt
		-t	--type		Input type of the Pangenome. Currently supports "Pangenome" <default>, Pangenome "Iterative", and "Panoct".

		
### make_conf_file.pl 
make_conf_file.pl generates a functional configuration file used by make_panacea.pl to
color the genes and make the tables for the terms. It is designed to work on the output
of the JCVI Pangenome pipeline.

**running make_conf_file.pl**
	
	make_conf_file.pl [options]
	Options:
        -d	--dir input pangenome directory. Default is current directory
		-o	--out output directory. Default is current directory
        -g	--go location of go term obo file. Default is current directory
        -a	--aro location of aro obo file. Default is current directory
        -h 	--help 


### make_rgi_clusters.pl 
make_rgi_clusters.pl translates the RGI output file dataSummary.txt into the PanGenome annotation
table format. Required to make it readable for PanACEA. Either best or all hits can be outputed

**running make_conf_file.pl**
	
	make_rgi_clusters.pl [options] -d directory -t type -o outputFlatFile
	Options:
        -i	--input Location of the dataSummary.txt file. Required
		-o	--output output file. Default is "aro_centroid.list.txt" in the current directory
        -t 	--type Output type of the file, either "best" (only the best hit) or "all". Best is default 

## Included Other Files

### PanACEA.manual.pdf
PDF manual including screen shots and usage guides for both the command line and the web
interface for PanACEA

### run_panacea_example.sh
Shell script to run the scripts required to build the PanACEA web pages for the examplar data

### example_dir/
Directory containing the example PanACEA data derived from the PanOCT run on the Enterobactor
hormaechei dataset from *Chavda et al 2016*. Can be generated using the command:

	run_panacea_example.sh <output_directory>

### example_dir/sample_output.ehormaechei.tar.gz
Compressed file containing the multi-file output of the PanACEA run on Enterobactor. 
