#!/usr/bin/python
from sys import argv
import os
import pandas as pd

# To be able to build a custom Kraken2 database, the fasta files need to have taxids added to the header.
# This script will read a fasta file, obtain the taxid from a tsv file,
# and add the taxid to the header of the fasta file as required by Kraken2. 
# It also generates an accession2taxid file that maps the contig accession numbers to the taxids, this file is needed for Kraken2 to work properly, and if genomes are selected to build a custom database, the file should include all the contig accessions in the database.

# Note that this script works on fastas downloaded from NCBI, which have the accession number in the name of the fasta file itself.

# The script can be used in bash as follows:
# taxids_added_dir="/path/to/output/"
# for file in *.fna
#    do python3 add_taxids_for_kraken2.py $file "$taxids_added_dir"$(echo $(basename $file .fna))_taxids.fna
# done


# Function to create the output fasta file path
def create_out_fasta_path(fasta_file):
    path_to_out_folder="/path/to/output/"
    base=os.path.basename(fasta_file)
    filename = os.path.splitext(base)[0]
    out_fasta_file = path_to_out_folder + filename + "_taxids_added.fna"
    return out_fasta_file

# Function to get the accession number from the name of the fasta file, and get the taxid from the tsv file
def obtain_taxid(fasta_file, acc_taxids_file):
    path = os.path.basename(fasta_file)
    # accession is the first and second part of the name of the fasta file if it is separated by _
    acc = "_".join(path.split("_")[:2])
    taxids_df = pd.read_csv(acc_taxids_file)
    # Create a dictionary with accession numbers as keys and taxids as values
    acc_taxid_dict = taxids_df.set_index(taxids_df.columns[0])[taxids_df.columns[2]].to_dict()
    taxid = acc_taxid_dict.get(acc)
    return taxid

# Function to parse the fasta and add the taxid as needed by Kraken2
def parse_fasta(fasta_file, taxid, out_fasta):
    with open(fasta_file) as f1:
        with open(out_fasta, "w") as fasta_taxids:
            for line in f1:
                if line.startswith(">"):
                    fasta_taxids.write(line.strip() + "|kraken:taxid|" + str(taxid) + "\n")
                else:
                    fasta_taxids.write(line.strip() + "\n")

# Function to get the list of contig accessions from the fasta file
def get_contig_accessions(fasta_file):
    chrom_accession_list=[]
    with open(fasta_file) as f1:
        for line in f1:
            if line.startswith(">"):
                chrom_accession_list.append(line.strip().split(" ")[0][1:])
    return (chrom_accession_list)

# Function to save the accession to taxid mapping in a file
def save_acc2taxid(chrom_accession_list,acc_taxids_file,taxid):
    with open(acc_taxids_file,"a") as out_acc2taxids:
        for i in chrom_accession_list:
            out_acc2taxids.write(i.split(".")[0] + "\t" + i + "\t" + str(taxid) + "\t" + "na" + "\n")

# Main function to run the script
def main():
    fasta_file = argv[1]
    out_fasta_file = create_out_fasta_path(fasta_file)
    acc_taxids_file = "/path/to/file/with/accessions_and_taxids.csv"
    taxid = obtain_taxid(fasta_file,acc_taxids_file)
    parse_fasta(fasta_file, taxid, out_fasta_file)
    acc_2_taxid = "/path/to/kraken2/database/DATABASENAME.accession2taxid"
    save_acc2taxid(get_contig_accessions(fasta_file),acc_2_taxid,taxid)

if __name__ == "__main__":
    main()
