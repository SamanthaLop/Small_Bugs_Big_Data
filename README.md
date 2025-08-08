# Small_Bugs_Big_Data
Selected custom scripts used in the analysis of Malaise trap bulkDNA shotgun sequencing data, as part of the manuscript "Small Bugs, Big Data: Metagenomics for arthropod biodiversity monitoring" - for submission to Ecology and Evolution.

______

### 1. add_taxids_for_kraken2.py
This is a python script that will add the corresponding taxIDs to each fasta header in a fasta file using the name of the fasta file (which is assumed to include the accession for that assembly) and a tsv file that contains three rows: Assembly Accession, Organism Name, and Organism Taxonomic ID, where the latter refers to the NCBI taxID. The script will add the taxid as needed by [Kraken2](https://github.com/DerrickWood/kraken2) (`kraken:taxid|123456`) and generate the information required at the contig level to update the `accession2taxid` file required by Kraken2 for building a database.


### 2. kraken2_report_parsing_functions.py
This python script (used as a module) contains the functions used to parse kraken2 reports for: confirming the input path contains Kraken2 reports, extracting - for all reports in the input directory - the number of classified and unclassified reads (and percentages), calculating the total number of reads, and concatenating all kraken2 reports into a single dataframe for further analysis.


### 3. mitochondrial_allele_frequencies.py
This is a custom script that contains the functions used to read and plot allele frequencies from a dp4 file. The dp4 file is assumed to be generated through bcftools mpileup (with flag -a DP4) on a bam file which is then piped to bcftools call, resulting in an output_mpileup file which is then indexed and gzipped. The resulting file is then queried with `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP4]\n' ${output_mpileup}.gz > dp4_OUTPUT.tsv`. The functions in this script parse the dp4 file into a dataframe that is then used to compute allele frequencies, save them as a csv, and plot them as a histogram.
