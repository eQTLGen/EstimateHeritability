#!/usr/bin/env python3

"""
Created:      15/04/2023
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2023 C.A. Warmerdam

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import os
import sys
import argparse
import gzip

import numpy as np
import pandas as pd

# Metadata
__program__ = "Annotate Loci"
__author__ = "C.A. (Robert) Warmerdam"
__email__ = "c.a.warmerdam@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


# Constants

# Classes
class MafCalculator:
    def __init__(self, inclusion_path, cohorts, maf_table, flipped, table_name="filter_logs.log",
                 variant_inclusion_format="%s_SnpsToInclude.txt", gene_inclusion_format="%s_GenesToInclude.txt"):
        self.cohorts = cohorts
        overview_df = pd.read_table(os.path.join(inclusion_path, table_name), index_col=False)
        overview_df.set_index('Dataset', inplace=True)
        print(overview_df.head())
        self.overview_df = overview_df.loc[self.cohorts]
        self.maf_table = maf_table[self.overview_df.index].copy()
        self.maf_table.loc[flipped, :] = 1 - self.maf_table.loc[flipped, :]
        self.overview_df['snp_inclusion_path'] = (
            self.overview_df.index.map(lambda name: os.path.join(inclusion_path, variant_inclusion_format % name)))
        self.overview_df['gene_inclusion_path'] = (
            self.overview_df.index.map(lambda name: os.path.join(inclusion_path, gene_inclusion_format % name)))
        self.snp_inclusion_df = self.load_inclusion_df('snp_inclusion_path')
        self.gene_inclusion_df = self.load_inclusion_df('gene_inclusion_path')
    def load_inclusion_df(self, column):
        # Create an empty dictionary to store dataframes
        dfs = {}
        # Generate dict of inclusion paths
        inclusion_paths = self.overview_df[column].to_dict()
        # Iterate over files in the directory
        for cohort, filepath in inclusion_paths.items():
            print(cohort, filepath)
            df = pd.read_csv(filepath)
            print(df.head())
            df[cohort] = 1
            # Set the index to the 'ID' column
            df.set_index('ID', inplace=True)
            # Add the dataframe to the dictionary
            dfs[cohort] = df
        # Merge the dataframes on their index
        merged_df = pd.concat(dfs.values(), axis=1)
        # Replace NaN values with 0 and convert to boolean
        merged_df = merged_df.fillna(0).astype(bool)
        return merged_df
    def calculate_maf(self, gene_variant_df):
        # Reformat the presence of the variants in the given dataframe
        variant_presence = (
            gene_variant_df
            .merge(self.snp_inclusion_df, right_index=True, left_on='variant', how='left')
            .set_index(['phenotype', 'variant']))
        print(variant_presence)
        # Reformat the presence of the genes in the given dataframe
        gene_presence = (
            gene_variant_df
            .merge(self.gene_inclusion_df, right_index=True, left_on='phenotype', how='left')
            .set_index(['phenotype', 'variant']))
        print(gene_presence)
        # Now that the tables displaying presence have both the same index, we can determine
        # for each combination if it is present or not.
        combined_presence = variant_presence & gene_presence
        print(combined_presence)
        # Now, reformat the maf table to also be according to this format.
        variant_maf = (
            gene_variant_df
            .merge(self.maf_table, right_index=True, left_on='variant', how='left')
            .set_index(['phenotype', 'variant']))
        # Now, for each cohort in the maf table, multiply all MAFs by the sample size
        variant_maf_weighted = variant_maf.mul(self.overview_df.loc[variant_maf.columns, "N"])
        print(variant_maf_weighted)
        # Now sum the weighted MAFs, and divide this by the total sample size
        maf = (variant_maf_weighted.where(combined_presence).sum(axis=1)
               / combined_presence.dot(self.overview_df["N"]))
        # Now return this
        return maf


class GencodeParser:
    def __init__(self, filename):
        self.filename = filename
        self.df = self.parse()
    def parse(self):
        genes = []
        with gzip.open(self.filename, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if fields[2] == 'gene':
                        attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                        genes.append({'gene_id': attributes['gene_id'],
                                      'chromosome': fields[0],
                                      'start': int(fields[3]),
                                      'end': int(fields[4])})
        df = pd.DataFrame(genes).astype({'start': 'Int64', 'end': 'Int64'})
        df['chromosome'] = df['chromosome'].str.lstrip('chr').replace({"M": "25", "X": "23", "Y": "24"}).astype('Int64')
        df['gene_id'] = df['gene_id'].str.split('.').str[0]
        print(df['chromosome'].unique())
        return df[['gene_id', 'chromosome', 'start', 'end']].drop_duplicates(keep='first', subset=["gene_id"])


class GtfParser:
    CHROMOSOMES = np.concatenate([
        np.array(list(range(1,26))).astype('str'),
        np.array(["M", "MT", "X", "Y"])])
    def __init__(self, filename):
        self.filename = filename
        self.df = self.parse()
    def parse(self):
        genes = []
        with gzip.open(self.filename, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if fields[2] == 'gene':
                        attributes = dict(item.strip().split(' ') for item in fields[8].split(';') if ' ' in item.strip())
                        genes.append({'gene_id': attributes['gene_id'].strip("\""),
                                      'chromosome': fields[0],
                                      'start': int(fields[3]),
                                      'end': int(fields[4])})
        df = pd.DataFrame(genes).astype({'start': 'Int64', 'end': 'Int64'})
        df = df.loc[df['chromosome'].isin(self.CHROMOSOMES)]
        df['chromosome'] = df['chromosome'].str.lstrip('chr')\
            .replace({"M": "25", "MT": "25", "X": "23", "Y": "24"}).astype('Int64')
        df['gene_id'] = df['gene_id'].str.split('.').str[0]
        print(df['chromosome'].unique())
        return df[['gene_id', 'chromosome', 'start', 'end']].drop_duplicates(keep='first', subset=["gene_id"])


# Functions

# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = argparse.ArgumentParser(description='Annotate significant eQTL results with variant and gene information')

    parser.add_argument('--input-file', dest='input_file', required=True,
                        help='Path to the table containing eQTL results')
    parser.add_argument('--variant-reference', dest='variant_reference', required=True,
                        help='Path to the table containing all SNPs from a reference panel')
    parser.add_argument('--gene-ref', dest='gene_ref', required=True,
                        help='Path to a GTF or GFF3 file')
    parser.add_argument('--out-prefix', dest='out_prefix', required=True,
                        help='Prefix to use for output file names')
    parser.add_argument('--cohorts', dest='cohorts', required=True, nargs='+',
                            help='Names of cohorts used in the meta-analysis')
    parser.add_argument('--inclusion-path', dest='inclusion_path', required=True,
                        help='Inclusion_path')

    args = parser.parse_args(argv[1:])
    print(args)
    print("Loading variant reference from '{}'".format(args.variant_reference))

    variant_reference = (
        pd.read_csv(args.variant_reference, sep=' ', dtype={'CHR': "Int64", 'bp': "Int64"})
        .drop(["allele1", "allele2"], axis=1)
        .rename({"ID": "variant", "bp": "bp_variant", "CHR": "chromosome_variant",
                 "str_allele1": "allele_ref", "str_allele2": "allele_eff"}, axis=1))

    print("Variant reference loaded:")
    print(variant_reference.head())

    print("Loading gene annotations from '{}'".format(args.gene_ref))

    gencode_parser = GtfParser(args.gene_ref)
    gene_dataframe = gencode_parser.df

    print(gene_dataframe.head())

    overview_df = pd.read_table(os.path.join(args.inclusion_path, "filter_logs.log"), index_col=False)
    overview_df.set_index('Dataset', inplace=True)
    total_sample_size = overview_df.loc[np.array(args.cohorts), "N"].sum()

    eqtls = pd.read_csv(args.input_file, sep='\t')

    print(eqtls.head())

    i_squared_threshold = 40
    sample_size_threshold = total_sample_size * 0.5
    pass_sample_size_threshold = (eqtls.sample_size > sample_size_threshold)
    pass_i_squared_threshold = eqtls.i_squared < i_squared_threshold

    passed_variants = np.logical_and(pass_sample_size_threshold, pass_i_squared_threshold)

    print("Total theoretical sample size = {}".format(total_sample_size))
    print("Sample size threshold = {}".format(sample_size_threshold))
    print("Frequency of passed variants:")
    print(np.unique(pass_sample_size_threshold, return_counts=True))

    print("I2 threshold = {}".format(i_squared_threshold))
    print("Frequency of passed variants:")
    print(np.unique(pass_i_squared_threshold, return_counts=True))

    print("Filtering variants")
    print("Frequency of passed variants:")
    print(np.unique(passed_variants, return_counts=True))

    eqtls_filtered = eqtls.loc[passed_variants]
    # For each gene, check what the number of variants are. If it is lower than

    if eqtls_filtered.shape[0] < 0.5 * eqtls.shape[0]:
        print("The number of filtered variants is under 50% of the number of input variants ({} vs {}, {}%)"
              .format(eqtls_filtered.shape[0], eqtls.shape[0], eqtls_filtered.shape[0]/eqtls.shape[0]))
        print("exiting...")
        return 0

    # Perform method
    eqtls_annotated = (
        eqtls_filtered
        .merge(variant_reference, how="inner", on="variant", validate="m:1")
        .merge(gene_dataframe, how="inner", left_on="phenotype", right_on="gene_id",
               suffixes=('', '_gene'), validate="m:1"))

    print(eqtls_annotated.head())

    # Identify genes that have a cis-effect
    cis = np.logical_and(
        eqtls_annotated.chromosome_variant == eqtls_annotated.chromosome_gene,
        np.logical_or((eqtls_annotated.bp_variant - eqtls_annotated.start).abs() < 1*10**6,
                      (eqtls_annotated.bp_variant - eqtls_annotated.end).abs() < 1*10**6))

    trans = ~np.logical_and(
        eqtls_annotated.chromosome_variant == eqtls_annotated.chromosome_gene,
        np.logical_or((eqtls_annotated.bp_variant - eqtls_annotated.start).abs() < 5*10**6,
                      (eqtls_annotated.bp_variant - eqtls_annotated.end).abs() < 5*10**6))

    ldsc_selector = {"variant": "SNP", "sample_size": "N", "z_score": "Z", "p_value": "P",
                     "beta": "BETA", "standard_error": "SE", "allele_eff": "A1", "allele_ref": "A2"}

    ldsc_processed_output = eqtls_annotated.rename(columns=ldsc_selector)[[*ldsc_selector.values()]]

    if ldsc_processed_output.loc[cis].shape[0] > 0:
        ldsc_processed_output.loc[cis].to_csv("{}_cis.csv.gz".format(args.out_prefix), sep="\t", index=False)
    if ldsc_processed_output.loc[trans].shape[0] > 0:
        ldsc_processed_output.loc[trans].to_csv("{}_trans.csv.gz".format(args.out_prefix), sep="\t", index=False)
    if ldsc_processed_output.shape[0] > 0:
        ldsc_processed_output.to_csv("{}_gw.csv.gz".format(args.out_prefix), sep="\t", index=False)

    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
