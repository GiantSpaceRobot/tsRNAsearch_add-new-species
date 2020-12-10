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
gunzip *.gz
# or
pigz -p 5 *.gz

# Run as follows:
Usage: $0 -s mouse -f mm10-tRNAs.fa -A Mus_musculus.GRCm38.95.gtf -F Mus_musculus.GRCm38.dna.primary_assembly.fa -o OutputDirectory
" 1>&2; }
info() { echo "
Options

	-h	Print the usage and options information
	-s	Species (e.g. mouse)
	-f	FASTA file of tRNAs for new species
	-A	Full GTF file for new species
	-F	Full FASTA file containing genome of new species 
	-o	Output directory for the new GTFs/FASTA files

	" 1>&2; }

while getopts ":hs:f:A:F:o:" o; do
    case "${o}" in
		h)
			usage
			info
			exit 0
			;;
		s)
			species="$OPTARG"
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

echo "Creating new directories..."
mkdir -p $outDir
mkdir -p $outDir/Intermediate-files

tRNA_basename="$(basename -- $in_tRNAs)"
tRNA_newname_FA="${tRNA_basename}_relative_temp.fa"
tRNA_newname_GTF="${tRNA_basename}_relative_temp.gtf"
ncRNA_basename="${species}_ncRNAs.gtf"
ncRNA_newname_FA="${ncRNA_basename}_relative.fa"
ncRNA_newname_GTF="${ncRNA_basename}_relative.gtf"
genome_basename="$(basename -- $in_genome)"

echo "Extracting features from GTF..."
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
	> $outDir/Intermediate-files/ncRNAs.gtf 
grep Mt_tRNA \
	$in_ncRNAs_GTF \
	> $outDir/Intermediate-files/Mt_tRNAs.gtf

### Get intron coordinates from .ss file
#grep -Ev ^'HMM|Seq:|Str:|    ' $in_tRNAs_ss \
#	| tr '\n' '\t' \
#	| sed 's/chr/\nchr/g' \
#	| grep intron \
#	| unexpand \
#	| awk '{print $1"-"$7$9"\t"$17}' \
#	| rev \
#	| sed 's/-/\t/' \
#	| rev \
#	| sort -V \
#	> $outDir/"${species}_tRNA-introns-for-removal.tsv"

### Create new tRNA and ncRNA GTFs and FASTAs
echo "Creating new GTF and FASTA files..."
python bin/FASTA-to-GTF.py tRNA $in_tRNAs $outDir/Intermediate-files/tRNA_newname.fa $outDir/Intermediate-files/tRNA_newname &    
python bin/GFF3-to-FASTA.py $outDir/Intermediate-files/ncRNAs.gtf $in_genome $outDir/Intermediate-files/ncRNAs_relative.fa $outDir/Intermediate-files/ncRNAs_relative.gtf &

wait

python bin/tRNA-GTF-to-FASTA.py $outDir/Intermediate-files/tRNA_newname_absolute.gtf $in_genome $outDir/Intermediate-files/tRNA_intermediate.fa $outDir/Intermediate-files/tRNA_intermediate_relative.gtf 

### Run cd-hit using similarity cut-off of 99.5%
echo "Running cd-hit on ncRNAs..."
cd-hit-est -i $outDir/Intermediate-files/ncRNAs_relative.fa -o $outDir/Intermediate-files/ncRNAs_relative_cdhit.fa -c 0.995
#cd-hit-est -i $tRNA_newname_FA -o "${tRNA_basename}_relative_cdhit.fa" -c 0.995

### Combine the tRNA and ncRNA FASTAs
echo "Combining ncRNA and tRNA FASTA files..."
cat $outDir/Intermediate-files/ncRNAs_relative_cdhit.fa $outDir/Intermediate-files/tRNA_intermediate.fa > $outDir/"${species}_tRNAs-and-ncRNAs_relative.fa" &

### Run Reduce-GTF - This removes multiple GTF lines referring to the same feature
echo "Removing GTF duplicates..."
python bin/Reduce-GTF.py $outDir/Intermediate-files/ncRNAs_relative_cdhit.fa $outDir/Intermediate-files/ncRNAs_relative.gtf $outDir/Intermediate-files/"${species}_ncRNAs_relative_cdhit.gtf" & 
python bin/Reduce-GTF.py $outDir/Intermediate-files/tRNA_intermediate.fa $outDir/Intermediate-files/tRNA_intermediate_relative.gtf $outDir/"${species}_tRNAs_relative.gtf" &

wait

sed 's/gene_type/gene_biotype/g' $outDir/Intermediate-files/"${species}_ncRNAs_relative_cdhit.gtf" > $outDir/"${species}_ncRNAs_relative_cdhit.gtf" & # Make sure ncRNA type is defined using gene_biotype

### Find tRNA introns for removal
#awk '{print $1}' $outDir/"${species}_tRNAs_relative_cdhit.gtf" | uniq -d > $outDir/Intermediate-files/"${species}_tRNAs-with-introns.txt"
#join <(sort $outDir/Intermediate-files/"${species}_tRNAs-with-introns.txt") <(sort $outDir/"${species}_tRNAs_relative_cdhit.gtf") > $outDir/Intermediate-files/"${species}_Intron-containing-tRNAs.gtf"
#python bin/Remove-intron.py $outDir/Intermediate-files/"${species}_Intron-containing-tRNAs.gtf" $outDir/"${species}_tRNA-introns-for-removal.tsv"

### Get names of all tRNAs and ncRNAs
echo "Gathering all ncRNA and tRNA names..."
grep '>' $outDir/"${species}_tRNAs-and-ncRNAs_relative.fa" | \
	sed 's/>//g' | \
	awk -F ' ' '{print $1}' \
	> $outDir/"${species}_all-ncRNAs.txt"


### Create empty count files
echo "Creating empty count files..."
sed 's/$/\t0/' $outDir/"${species}_all-ncRNAs.txt" > $outDir/Intermediate-files/"${species}_all-ncRNAs.count" # Add tab and 0 after every ncRNA
grep -i tRNA $outDir/Intermediate-files/"${species}_all-ncRNAs.count" > $outDir/"${species}_empty_tRNA.count" # Get tRNA count file
grep -iv tRNA $outDir/Intermediate-files/"${species}_all-ncRNAs.count" > $outDir/"${species}_empty_ncRNA.count" # Get ncRNA count file

### Get lengths of all tRNAs
echo "Generating numbers of and lengths of tRNAs..."
python bin/tRNA-lengths.py $outDir/"${species}_tRNAs_relative.gtf" $outDir/"${species}_tRNA-lengths.txt"

echo "ncRNAs in input: "$(grep -c '>' $outDir/Intermediate-files/ncRNAs_relative.fa)
echo "ncRNAs in output: "$(grep -c '>' $outDir/Intermediate-files/ncRNAs_relative_cdhit.fa)
echo "tRNAs in input: "$(grep -c '>' $in_tRNAs)
echo "tRNAs in output: "$(grep -c '>' $outDir/Intermediate-files/tRNA_intermediate.fa)

echo "Finished"

