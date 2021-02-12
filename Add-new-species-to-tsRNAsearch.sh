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
wget -q http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/tRNAs.fa .

# Gunzip files
gunzip *.gz
# or
pigz -p 5 *.gz

# Run as follows:
Usage: $0 -s ${species} -f tRNAs.fa -A Mus_musculus.GRCm38.95.gtf -F Mus_musculus.GRCm38.dna.primary_assembly.fa -o OutputDirectory
" 1>&2; }
info() { echo "
Options

	-h	Print the usage and options information
	-s	Species (e.g. ${species})
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
cat \
	$outDir/Intermediate-files/ncRNAs.gtf \
	$outDir/Intermediate-files/Mt_tRNAs.gtf \
	> $outDir/Intermediate-files/All-ncRNA_absolute-coordinates.gtf

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

### Add up- and downstream nucleotides to tRNAs
# Nuclear tRNAs
python bin/tRNA-GTF-to-FASTA.py $outDir/Intermediate-files/tRNA_newname_absolute.gtf $in_genome $outDir/Intermediate-files/tRNA_intermediate.fa $outDir/Intermediate-files/tRNA_intermediate_relative.gtf 
# Mitochondrial tRNAs
python bin/tRNA-GTF-to-FASTA.py $outDir/Intermediate-files/Mt_tRNAs.gtf $in_genome $outDir/Intermediate-files/Mt-tRNA_intermediate.fa $outDir/Intermediate-files/Mt-tRNA_intermediate_relative.gtf 

### Combine tRNA files
cat $outDir/Intermediate-files/tRNA_intermediate.fa $outDir/Intermediate-files/Mt-tRNA_intermediate.fa \
	> $outDir/Intermediate-files/All-tRNAs_intermediate.fa
cat $outDir/Intermediate-files/tRNA_intermediate_relative.gtf $outDir/Intermediate-files/Mt-tRNA_intermediate_relative.gtf \
	> $outDir/Intermediate-files/All-tRNAs_intermediate_relative.gtf

### Run cd-hit using similarity cut-off of 99.5%
echo "Running cd-hit on ncRNAs..."
cd-hit-est -i $outDir/Intermediate-files/ncRNAs_relative.fa -o $outDir/Intermediate-files/ncRNAs_relative_cdhit.fa -c 0.995
#cd-hit-est -i $tRNA_newname_FA -o "${tRNA_basename}_relative_cdhit.fa" -c 0.995

### Combine the tRNA and ncRNA FASTAs
echo "Combining ncRNA and tRNA FASTA files..."
cat $outDir/Intermediate-files/ncRNAs_relative_cdhit.fa $outDir/Intermediate-files/All-tRNAs_intermediate.fa > $outDir/"${species}_tRNAs-and-ncRNAs_relative.fa" &

### Run Reduce-GTF - This removes multiple GTF lines referring to the same feature
echo "Removing GTF duplicates..."
python bin/Reduce-GTF-v2.py $outDir/Intermediate-files/ncRNAs_relative_cdhit.fa $outDir/Intermediate-files/ncRNAs_relative.gtf $outDir/Intermediate-files/"${species}_ncRNAs_relative_cdhit.gtf" & 
python bin/Reduce-GTF-v2.py $outDir/Intermediate-files/All-tRNAs_intermediate.fa $outDir/Intermediate-files/All-tRNAs_intermediate_relative.gtf $outDir/"${species}_tRNAs_relative.gtf" &

wait

### Change "gene_type" to "gene_biotype" for consistency 
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
echo "tRNAs (including mitochondrial tRNAs) in output: "$(grep -c '>' $outDir/Intermediate-files/tRNA_intermediate.fa)


##### Add lookalikes to final FASTA file and GTF #####
echo "=====Adding ncRNA/tRNA lookalikes to FASTA file====="
### Compare tRNAs to genome
echo "Comparing tRNAs to genome..."
mmseqs easy-search --min-aln-len 15 --search-type 3 ${in_tRNAs} ${in_genome} $outDir/Intermediate-files/mmseqs_tRNAs_vs_${species}-genome.m8 tmpdir

### Get top hits for above comparison (alignment over 65 nt and identity = 100%)
echo "Getting top hits..."
awk '$3 == "1.000" {print $0}' $outDir/Intermediate-files/mmseqs_tRNAs_vs_${species}-genome.m8 | awk '$4 > 65 {print $0}' | sort -k1,1n > $outDir/Intermediate-files/mmseqs_tRNAs_vs_${species}-genome_top-hits.m8 

### Pull chromosome, start and stop coordinates from the top tRNA-vs-genome hits
echo "Extracting top hit coordinates..."
awk '{print $2"\t"$9"\t"$10}' $outDir/Intermediate-files/mmseqs_tRNAs_vs_${species}-genome_top-hits.m8 | sort -k1,1n | uniq > $outDir/Intermediate-files/Chrom_start_stop.tsv 

### Create genome FASTA masking tRNA top hits and all ncRNAs from GTF (i.e. replace all ncRNAs with NNNs in genome)
echo "Replace tRNAs/ncRNAs in genome with ambiguous nucleotides (NNNs)"
python bin/GTF-to-FASTA_EverythingButGTF-v3.py $outDir/Intermediate-files/Chrom_start_stop.tsv $outDir/Intermediate-files/All-ncRNA_absolute-coordinates.gtf $in_genome $outDir/Intermediate-files/${species}_genome_no-ncRNAs.fa 

### Compare tRNAs + ncRNAs against masked genome (i.e. genome excluding all perfect matches to ncRNAs)
### The aim here is to generate a list of ncRNA lookalike
echo "Compare ncRNAs/tRNAs to masked genome..."
mmseqs easy-search --search-type 3 $outDir/${species}_tRNAs-and-ncRNAs_relative.fa $outDir/Intermediate-files/${species}_genome_no-ncRNAs.fa $outDir/Intermediate-files/Lookalikes.m8 tmp 

### Convert lookalikes to FASTA sequence and create accompanying GTF
echo "Generate FASTA sequences of lookalikes..."
python2 bin/M8-to-FASTA.py $outDir/Intermediate-files/Lookalikes.m8 $outDir/Intermediate-files/${species}_genome_no-ncRNAs.fa $outDir/Intermediate-files/Lookalikes.fa $outDir/Intermediate-files/Lookalikes.gtf 

### Remove redundant sequences in lookalikes FASTA
echo "Remove redundant sequences from lookalikes..."
cd-hit-est -i $outDir/Intermediate-files/Lookalikes.fa -o $outDir/Intermediate-files/Lookalikes_cdhit.fa -c 0.995 

### Reduce GTF to match new reduced FASTA
echo "Reduce lookalike GTF to match reduced lookalike FASTA..."
python bin/Reduce-GTF-v2.py $outDir/Intermediate-files/Lookalikes_cdhit.fa $outDir/Intermediate-files/Lookalikes.gtf $outDir/Intermediate-files/Lookalikes_cdhit.gtf 

### Combine lookalikes with regular FASTA
echo "Combined ncRNA/tRNA lookalikes with genuine ncRNA/tRNA FASTA"
cat $outDir/"${species}_tRNAs-and-ncRNAs_relative.fa" \
	$outDir/Intermediate-files/Lookalikes_cdhit.fa \
	> $outDir/"${species}_tRNAs-and-ncRNAs-and-lookalikes_relative.fa"

echo "Finished"
