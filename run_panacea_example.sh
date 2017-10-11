#!/bin/sh
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


SCRIPT=$(readlink -f $0)
SRC_DIR=$(dirname $SCRIPT)
PAN_DIR="$SRC_DIR/example_dir/"
PAN_DIR="$(cd $PAN_DIR; pwd)"
TMP_DIR=$PAN_DIR
if [[ $# -eq 1 ]]; then
	TMP_DIR=$1	
fi
if [ ! -d "$TMP_DIR" ]; then 
	mkdir $TMP_DIR
fi

if [ -e "$PAN_DIR/aro_centroid.list.txt" ];
then
	echo "Found ARO mapfile"
else
	echo "Cannot find ARO mapfile. Making it..."
	perl $SRC_DIR/make_rgi_clusters.pl -i $PAN_DIR/dataSummary.txt -o $PAN_DIR/aro_centroid.list.txt
fi

TMP_DIR="$(cd $TMP_DIR; pwd)"
echo "Making PanACEA function file configure file ..."
perl $SRC_DIR/make_conf_file.pl -d $PAN_DIR -o $TMP_DIR
TREE_FILE="$PAN_DIR/genome.tree"
FASTA_FILE="$PAN_DIR/combined.fasta"
echo "Making PanACEA HTML files ..."
perl $SRC_DIR/make_panacea.pl -i $PAN_DIR/PanACEA.flatfile.txt -o $TMP_DIR/ -f $TMP_DIR/func_file.conf.txt -t $TREE_FILE -a $FASTA_FILE -m $PAN_DIR/panacea_example.grouping
echo "Finished.."
