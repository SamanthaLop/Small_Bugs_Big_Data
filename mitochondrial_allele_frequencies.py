#!/usr/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# This script contains the functions to read and plot mitochondrial allele frequencies from a dp4 file.
# The DP4 file is generated with bcftools mpileup and flag -a DP4 on a bam file which is the piped to bcftools call into an OUTPUT_MPILEUP file that can be indexed and gzipped.
# The following bash command can be used to generate the DP4 file that is used in this script:
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP4]\n' ${output_mpileup}.gz > dp4_OUTPUT.tsv

# Function to read the DP4 file and return a DataFrame
def read_csv(file_dir):
    try:
        df = pd.read_csv(file_dir, sep="\t", header=None, names=["CHROM", "POS", "REF", "ALT", "DP4"])
        return df
    except Exception as e:
        print(f"Error reading {file_dir}: {e}")
        return None

# Function to split the DP4 column into four separate columns: REF_fwd, REF_rev, ALT_fwd, ALT_rev
def split_dp4_column(df):
    """
    Splits the DP4 column into four separate columns: REF_fwd, REF_rev, ALT_fwd, ALT_rev.
    """
    dp4_split = df["DP4"].str.split(",", expand=True).astype(int)
    df["REF_fwd"] = dp4_split[0]
    df["REF_rev"] = dp4_split[1]
    df["ALT_fwd"] = dp4_split[2]
    df["ALT_rev"] = dp4_split[3]
    return df

# Function to compute allele frequencies
def compute_allele_frequencies(df):
    """
    Computes allele frequencies from the DP4 columns.
    """
    df["REF_count"] = df["REF_fwd"] + df["REF_rev"]
    df["ALT_count"] = df["ALT_fwd"] + df["ALT_rev"]
    df["TOTAL_depth"] = df["REF_count"] + df["ALT_count"]
    # Avoid division by zero
    df["ALT_freq"] = df.apply(lambda row: row["ALT_count"] / row["TOTAL_depth"] if row["TOTAL_depth"] > 0 else 0, axis=1)
    return df

# Function to save the DataFrame to a CSV file
def save_to_csv(df, sample_name):
    """
    Saves the DataFrame to a CSV file.
    """
    output_file = f"{directory}/allele_freq_{sample_name}.tsv"
    df[["CHROM", "POS", "REF", "ALT", "REF_count", "ALT_count", "TOTAL_depth", "ALT_freq"]].to_csv(output_file, sep="\t", index=False)
    print(f"Allele frequencies saved to {output_file}")

# Function to plot allele frequencies
def plot_allele_frequencies(df, sample_name):
    """
    Plots the allele frequencies and saves the plot.
    """
    plt.figure(figsize=(10, 6))
    df["ALT_freq"].hist(bins=50, edgecolor='black')
    plt.xlabel("ALT Allele Frequency")
    plt.ylabel("Number of sites")
    plt.title(f"Histogram of ALT Allele Frequencies for {sample_name}")
    plt.grid(True)
    plt.tight_layout()
    # Save as SVG for later editing
    plt.savefig(f"{directory}/allele_freq_histogram_{sample_name}_5.svg", format='svg')
    # Save as PNG
    plt.savefig(f"{directory}/allele_freq_histogram_{sample_name}_5.png")
    plt.clf()
    plt.cla()
    plt.close()


# To run the functions, first state the directory, file_dir, and sample_name
directory = "path/to/your/directory/with/dp4/file"
file_dir = f"{directory}/dp4_OUTPUT.tsv"
sample_name = "YOUR_SAMPLE_NAME"

# Then you can use the functions above
df = read_csv(file_dir)
if df is not None:
    df = split_dp4_column(df)
    df = compute_allele_frequencies(df)
    save_to_csv(df, sample_name)
    plot_allele_frequencies(df, sample_name)