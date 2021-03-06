# Add new species to tsRNAsearch 

## Quickstart
All you need to do is: 
```
# Download the genome and GTF of the species you are interested in, and a tRNA FASTA file:
wget -q http://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz .
wget -q http://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz .
wget -q http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.fa .

# Gunzip files
gunzip *gtf.gz

# Run as follows (replace the word mouse with your species):
Usage: ./Add-new-species-to-tsRNAsearch.sh -s mouse -T mm10-tRNAs.fa -S Mus_musculus.GRCm38.95.gtf -G Mus_musculus.GRCm38.dna.primary_assembly.fa -o OutputDirectory
```
Move files to tsRNAsearch directory:
```
mv OutputDirectory/*.gtf /path/to/tsRNAsearch/DBs/
mv OutputDirectory/*tRNA-lengths.txt /path/to/tsRNAsearch/additional-files/
mv OutputDirectory/*all-ncRNAs.txt /path/to/tsRNAsearch/additional-files/
mv OutputDirectory/*tRNAs-and-ncRNAs_relative_cdhit.fa /path/to/tsRNAsearch/DBs/
mv OutputDirectory/*.count /path/to/tsRNAsearch/additional-files/
```
Add new species to tsRNAsearch (replace the word *mouse* with your species:
```
mkdir -p DBs/species_index/mouse-ncRNAs
STAR --runThreadN 1 --runMode genomeGenerate --genomeDir DBs/species_index/mouse-ncRNAs/ --genomeFastaFiles DBs/mouse_tRNAs-and-ncRNAs_relative_cdhit.fa
```

## More Information
Options:

* -h *Print the usage and options information*
* -s *Species (e.g. mouse)*
* -T *FASTA file of tRNAs for new species*
* -S *Full GTF file for new species*
* -G *Full FASTA file containing genome of new species* 
* -o *Output directory for the new GTFs/FASTA files*

## Contributors
* Paul Donovan, PhD (pauldonovandonegal@gmail.com)

## License
This project is licensed under the MIT License.

