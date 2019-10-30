#!/bin/bash
# Author: Paul Donovan 
# Email: pauldonovan@rcsi.com
# 23-10-19

usage() { echo "
==========
How to use
==========

# Download the genome and GTF of the species you are interested in, and a tRNA FASTA file:
wget -q http://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz .
wget -q http://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz .
wget -q http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.fa .

# Gunzip files
gunzip *gtf.gz

# Run as follows:
Usage: $0 -s mouse -a mm10-tRNAs-confidence-set.out -f mm10-tRNAs.fa -S Mus_musculus.GRCm38.95.gtf -G Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -o OutputDirectory
" 1>&2; }
info() { echo "
Options

	-h	Print the usage and options information
	-s	Species (e.g. mouse)
	-a	tRNA 'condfidence-set.out' file from GtRNAdb/tRNAscan
	-f	FASTA file of tRNAs for new species
	-A	Full GTF file for new species
	-F	Full FASTA file containing genome of new species 
	-o	Output directory for the new GTFs/FASTA files

	" 1>&2; }

while getopts ":hs:a:f:A:F:o:" o; do
    case "${o}" in
		h)
			usage
			info
			exit 0
			;;
		s)
			species="$OPTARG"
			;;
		a)
			in_tRNAs_ss="$OPTARG"
			;;
		f)
			in_tRNAs="$OPTARG"
			;;
		A)
			in_ncRNAs_GTF="$OPTARG"
			;;
		F)
			in_genome="$OPTARG"
			;;
		o)
			outDir="$OPTARG"
			;;
		*)
            echo "Error: incorrect input parameters!"
			usage
			info
			exit 1
            ;;
    esac
done

### If no command line arguments provided, quit
if [ -z "$*" ] ; then
	usage
    info
    echo "Error: no command line parameters provided!"
	exit 1
fi

mkdir -p $outDir
mkdir -p $outDir/Intermediate-files

tRNA_basename="$(basename -- $in_tRNAs)"
tRNA_newname_FA="${tRNA_basename}_relative.fa"
tRNA_newname_GTF="${tRNA_basename}_relative.gtf"
ncRNA_basename="${species}_ncRNAs.gtf"
ncRNA_newname_FA="${ncRNA_basename}_relative.fa"
ncRNA_newname_GTF="${ncRNA_basename}_relative.gtf"
genome_basename="$(basename -- $in_genome)"

### Get ncRNAs from GTF file
grep \
	-e miRNA \
	-e misc_RNA \
	-e rRNA \
	-e scRNA \
	-e snRNA \
	-e snoRNA \
	-e ribozyme \
	-e sRNA \
	-e scaRNA \
	$in_ncRNAs_GTF \
	> $outDir/Intermediate-files/$ncRNA_basename 

### Get intron coordinates from .ss file
grep -Ev ^'HMM|Seq:|Str:|    ' $in_tRNAs_ss \
	| tr '\n' '\t' \
	| sed 's/chr/\nchr/g' \
	| grep intron \
	| unexpand \
	| awk '{print $1"-"$7$9"\t"$17}' \
	| rev \
	| sed 's/-/\t/' \
	| rev \
	| sort -V \
	> $outDir/"${species}_tRNA-introns-for-removal.tsv"

### Create new tRNA and ncRNA GTFs and FASTAs
python bin/FASTA-to-GTF.py tRNA $in_tRNAs $tRNA_newname_FA $tRNA_newname_GTF &    
python bin/GFF3-to-FASTA.py $outDir/Intermediate-files/${ncRNA_basename} $in_genome $ncRNA_newname_FA $ncRNA_newname_GTF &

wait

### Run cd-hit using similarity cut-off of 99.5%
cd-hit-est -i $ncRNA_newname_FA -o "${ncRNA_basename}_relative_cdhit.fa" -c 0.995
cd-hit-est -i $tRNA_newname_FA -o "${tRNA_basename}_relative_cdhit.fa" -c 0.995

### Combine the tRNA and ncRNA FASTAs
cat "${ncRNA_basename}_relative_cdhit.fa" "${tRNA_basename}_relative_cdhit.fa" > $outDir/"${species}_tRNAs-and-ncRNAs_relative_cdhit.fa" &

### Run Reduce-GTF
python bin/Reduce-GTF.py "${ncRNA_basename}_relative_cdhit.fa" $ncRNA_newname_GTF $outDir/Intermediate-files/"${species}_ncRNAs_relative_cdhit.gtf" &
sed 's/gene_type/gene_biotype/g' $outDir/Intermediate-files/"${species}_ncRNAs_relative_cdhit.gtf" > $outDir/"${species}_ncRNAs_relative_cdhit.gtf" # Make sure ncRNA type is defined using gene_biotype
python bin/Reduce-GTF.py "${tRNA_basename}_relative_cdhit.fa" $tRNA_newname_GTF $outDir/"${species}_tRNAs_relative_cdhit.gtf" &

wait

mv "${tRNA_basename}_"* $outDir/Intermediate-files/
mv "${ncRNA_basename}_"* $outDir/Intermediate-files/

### Find tRNA introns for removal
#awk '{print $1}' $outDir/"${species}_tRNAs_relative_cdhit.gtf" | uniq -d > $outDir/Intermediate-files/"${species}_tRNAs-with-introns.txt"
#join <(sort $outDir/Intermediate-files/"${species}_tRNAs-with-introns.txt") <(sort $outDir/"${species}_tRNAs_relative_cdhit.gtf") > $outDir/Intermediate-files/"${species}_Intron-containing-tRNAs.gtf"
#python bin/Remove-intron.py $outDir/Intermediate-files/"${species}_Intron-containing-tRNAs.gtf" $outDir/"${species}_tRNA-introns-for-removal.tsv"

### Get names of all tRNAs and ncRNAs
grep '>' $outDir/"${species}_tRNAs-and-ncRNAs_relative_cdhit.fa" | \
	sed 's/>//g' | \
	awk -F ' ' '{print $1}' \
	> $outDir/"${species}_all-ncRNAs.txt"

### Create empty count files
sed 's/$/\t0/' $outDir/"${species}_all-ncRNAs.txt" > $outDir/Intermediate-files/"${species}_all-ncRNAs.count" # Add tab and 0 after every ncRNA
grep -i tRNA $outDir/Intermediate-files/"${species}_all-ncRNAs.count" > $outDir/"${species}_empty_tRNA.count" # Get tRNA count file
grep -iv tRNA $outDir/Intermediate-files/"${species}_all-ncRNAs.count" > $outDir/"${species}_empty_ncRNA.count" # Get ncRNA count file

### Get lengths of all tRNAs
python bin/tRNA-lengths.py $outDir/"${species}_tRNAs_relative_cdhit.gtf" $outDir/"${species}_tRNA-lengths.txt"

echo "Finished"

