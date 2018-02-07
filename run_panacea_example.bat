echo off
set dirID=%~dp0
perl %dirID%\make_panacea_flatfile.pl -i %dirID%\example_dir\ -o %dirID%\example_dir\PanACEA.flatfile.txt -t Pangenome
perl %dirID%\make_rgi_clusters.pl -i %dirID%\example_dir\dataSummary.txt -o %dirID%\example_dir/aro_centroid.list.txt
perl %dirID%\make_conf_file.pl -d %dirID%\example_dir\ -o %dirID%\example_dir\ -a %dirID%\example_dir\aro.obo -g %dirID%\example_dir\gene_ontology.1_2.obo
perl %dirID%\make_panacea.pl -i %dirID%\example_dir\PanACEA.flatfile.txt -o %dirID%\example_dir\ -f %dirID%\example_dir\func_file.conf.txt -t %dirID%\example_dir\genomes.tree -a %dirID%\example_dir\cluster_alignments,%dirID%\example_dir\combined.fasta