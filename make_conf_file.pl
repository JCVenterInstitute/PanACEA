#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use Cwd 'cwd';
use warnings;

my ($out_dir, $in_dir, $help, $aro_loc, $go_loc);
GetOptions("dir|d=s"=>\$in_dir, "out|o=s"=>\$out_dir, "help|h|?"=>\$help, "aro|a=s"=>\$aro_loc, "go|g=s"=>\$go_loc);
if ($help || !($in_dir) || !(-d $in_dir))
{
	print STDERR "Runs MASH of a directory or fasta file against the MASH sketch\n\n";
	print STDERR "--------------------USAGE--------------------------------------\n";
	print STDERR "	-d input pangenome directory\n";
	print STDERR "	-o output directory\n";
	print STDERR "	-g location of go term obo file. Default is input directory\n";
	print STDERR "	-a location of aro obo file. Default is input directory\n";
	print STDERR "	-h help\n";
		
	exit(0);
}
my @cwd = split "\/", abs_path($0);

my $src_dir = ""; for (my $i = 0; $i < scalar(@cwd)-1; $i++) { $src_dir .= $cwd[$i] . "/"; }

$in_dir = abs_path($in_dir);


my $num =2;
my @obo_files = ("DIR/aro.obo", "DIR/gene_ontology.1_2.obo");
my @types = ("AROTerms", "GOTerms");
my @map_files = ("DIR/aro_centroid.list.txt", "DIR/results/centroids.mapterm.txt");
my @color = ("ff0000", "CURR/cluster_colors.txt");
my @names = ("Antibiotic Resistance", "DIR/results/centroids.cluster_roles.txt");

$names[1] =~ s/DIR/$in_dir/g;
if (!-e $names[1])
{
	$names[1] = "DIR/centroids.cluster_roles.txt";
}
$map_files[1] =~ s/DIR/$in_dir/g;
if (!-e $map_files[1])
{
	$map_files[1] = "DIR/centroids.mapterm.txt";
	print STDERR "Here\n";

}

if ($aro_loc)
{
	if (-e $aro_loc)
	{
		$obo_files[0] = $aro_loc;
	}
	else
	{
		print STDERR "Cannot find aro obo file. Defaulting to input directory...";
	}
}

if ($go_loc)
{
	if (-e $go_loc)
	{
		$obo_files[1] = $go_loc;
	}
	else
	{
		print STDERR "Cannot find aro obo file. Defaulting to input directory...";
	}
}

if (!$out_dir)
{
	$out_dir = cwd();
}
if (!-e $out_dir)
{
	print STDERR "Cannot find directory. Quiting...\n\n\n"; exit();
}

if (-e "$out_dir/func_file.conf.txt")
{
	print STDERR "Functional file already exists... overwriting...\n";
}
open(OUT, ">", "$out_dir/func_file.conf.txt"); 

for (my $i = 0; $i < $num; $i++)
{
	if ($obo_files[$i] && $types[$i] && $map_files[$i] && $color[$i] && $names[$i] )
	{
		printf OUT "START\n";
		$types[$i] =~ s/DIR/$in_dir/g;
		$types[$i] =~ s/CURR/$src_dir/g;

		printf(OUT "ID\t%s\n",$types[$i]);
		$map_files[$i] =~ s/DIR/$in_dir/g;
		$map_files[$i] =~ s/CURR/$src_dir/g;

		printf(OUT "mapfile\t%s\n",$map_files[$i]);
		$obo_files[$i] =~ s/DIR/$in_dir/g;
		$obo_files[$i] =~ s/CURR/$src_dir/g;
		
		printf(OUT "ontology\t%s\n",$obo_files[$i]);
		$names[$i] =~ s/DIR/$in_dir/g;
		$names[$i] =~ s/CURR/$src_dir/g;

		printf(OUT "name\t%s\n",$names[$i]);
		$color[$i] =~ s/DIR/$in_dir/g;
		$color[$i] =~ s/CURR/$src_dir/g;

		printf(OUT "color\t%s\n",$color[$i]);
		printf(OUT "END\n\n")
	}
}
