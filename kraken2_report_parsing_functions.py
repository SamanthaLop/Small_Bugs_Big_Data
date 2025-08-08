#!/usr/bin/python
import pandas as pd
import sys
import os

# Function to check if the input directory is not empty and contains kraken2 reports
def check_k2_report_path(input_dir):
    """Checks if the input directory is not empty and contains kraken2 reports"""
    files=os.listdir(input_dir)
    assert (len(files)!=0), "The directory is empty"
    path=[]
    for f in files:
        fipath=os.path.join(input_dir,f)
        path.append(fipath)
        report=open(fipath, 'r')
        for line in report:
            fields=line.split('\t')
            cols=len(fields)
            assert (cols==8), "The %s file is not kraken2 summary report" %report
            break
    print ("Checking input file done")	
    return (path)

# Function to extract the number of classified and unclassified reads from kraken2 reports
def class_unclass_reads(input_dir):
    """Exracts, for each report, the number of classified and unclassified reads for each sample.
    The resulting dictionary has the following structure:
    (sample: [number_unclassified_reads, %_unclassified_reads, number_classified_reads, %_classified_reads])"""
    sample_dict={}
    for i in check_k2_report_path(input_dir):
        sample_dict[i] = None
        with open(i, "r") as report:
            h = []
            for line in report:
                fields = line.split("\t")
                if fields[5] == "U":
                    h.append(int(fields[1]))
                    h.append(float(fields[0].strip()))
                elif fields[5] == "R":
                    h.append(int(fields[1]))
                    h.append(float(fields[0].strip()))
        sample_dict[i]= h
    return sample_dict


# Function that will calculate the number of total reads by adding the number of classified reads and the number of unclassified reads from sample_dict
def calculate_total_reads(input_dict, sample_ids):
    """Calculates the number of total reads by adding the number of classified reads and the number 
    of unclassified reads from sample_dict, the output of class_unclass_reads
    """
    temp_dict = {}
    for i in input_dict:
        for j in sample_ids:
            if j in i:
                temp_dict[j] = input_dict[i][2] + input_dict[i][0]
    new_df = pd.DataFrame.from_dict(temp_dict, orient='index', columns=['Total'])
    return new_df


# Function that will make a dataframe using three columns from a dictionary:
# 1. sample_id, 2. number of classified reads, 3. number of unclassified reads
def number_reads_df(input_dict, sample_ids, percent):
    """Creates a dataframe with the sample id, the number of classified reads,
    and the number of unclassified reads. If percent == True, it will output
    values in percentage, if percent == False, it will output the true values
    """
    temp_dict = {}
    for i in input_dict:
        for j in sample_ids:
            if j in i and percent == False:
                temp_dict[j] = [input_dict[i][2], input_dict[i][0]]
            if j in i and percent == True:
                temp_dict[j] = [input_dict[i][3], input_dict[i][1]]
    new_df = pd.DataFrame.from_dict(temp_dict, orient='index', columns=['Classified', 'Unclassified'])
    return new_df

#  Function to get sample_ids
def get_sample_ids(input_tsv):
    """Creates a list of the sample ids
    """
    sample_ids = []
    for index,row in pd.read_csv(input_tsv, header=None).iterrows():
        sample_ids.append(row[0])
    return sample_ids

# Function to concatenate all kraken reports into one
def concatenate_kraken_files(input_dir, sample_dict):
    """Concatenates all kraken reports into one
    """
    all_data = []
    for file in sample_dict:
        data = pd.read_csv(file, sep='\t', header=None,
            names=["Percentage", "ReadsClade", "ReadsDirect",
            "NumMinimizers", "NumUniqueKmers", "Status",
            "TaxID", "Taxonomy"])
        data['Taxonomy'] = data['Taxonomy'].str.strip() 
        data.insert(0, "Sample", file.split("/")[-1].split(".txt")[0].split("_G")[0])
        all_data.append(data)
    concatenated_df = pd.concat(all_data)
    return concatenated_df

def calculate_total_reads(input_dict):
    temp_dict = {}
    for i in input_dict:
        for j in sample_ids:
            if j in i:
                temp_dict[j] = input_dict[i][2] + input_dict[i][0]
    new_df = pd.DataFrame.from_dict(temp_dict, orient='index', columns=['Total'])
    return new_df

# Function that will make a dataframe using three columns from a dictionary:
# 1. sample_id, 2. number of classified reads, 3. number of unclassified reads
def create_df_from_dict(input_dict, sample_ids, percent):
    temp_dict = {}
    for i in input_dict:
        for j in sample_ids:
            if j in i and percent == False:
                temp_dict[j] = [input_dict[i][2], input_dict[i][0]]
            if j in i and percent == True:
                temp_dict[j] = [input_dict[i][3], input_dict[i][1]]
    new_df = pd.DataFrame.from_dict(temp_dict, orient='index', columns=['Classified', 'Unclassified'])
    return new_df