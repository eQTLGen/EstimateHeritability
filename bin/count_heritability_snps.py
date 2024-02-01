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
__program__ = "Gene BED Files"
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
def count_heritability_snps_including(df, annot_columns, annot_file_df, frq_file_df):
    partial_inclusion = dict()
    for gene, grouped in df.groupby("name"):
        print(gene)
        sys.stdout.flush()
        # Get the rows from merged_annot_frq that do overlap with
        snp_inclusion = pd.concat([frq_file_df.query(
            "CHR == %s & BP > %d & BP < %d" %
            (row.chromosome, row.start, row.end))
            for index, row in grouped.iterrows()]).drop_duplicates()
        # Remove all duplicated rows and sum annotations for the merged rows,
        partial_inclusion[gene] = annot_file_df.iloc[
            snp_inclusion.index, annot_columns].sum()
    return pd.DataFrame(partial_inclusion)


def count_heritability_snps_excluding(df, annot_columns, annot_file_df, frq_file_df):
    partial = dict()
    for gene, grouped in df.groupby("name"):
        print(gene)
        sys.stdout.flush()
        # Get the rows from merged_annot_frq that do overlap with
        snp_inclusion = pd.concat([frq_file_df.query(
            "CHR == %s & BP > %d & BP < %d" %
            (row.chromosome, row.start, row.end))
            for index, row in grouped.iterrows()]).drop_duplicates()
        # Remove all duplicated rows and sum annotations for the merged rows,
        partial[gene] = annot_file_df.iloc[
            ~annot_file_df.index.isin(snp_inclusion.index), annot_columns].sum()
    return pd.DataFrame(partial)


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = argparse.ArgumentParser(
        description='Annotate significant eQTL results with variant and gene information')

    parser.add_argument('--bed-file-incl', dest='bed_inclusion_file', required=True, nargs="*",
                        help='Path to a .bed file')
    parser.add_argument('--bed-file-excl', dest='bed_exclusion_file', required=True, nargs="*",
                        help='Path to a .bed file')
    parser.add_argument('--out-prefix', dest='out_prefix', required=True,
                        help='Prefix to use for output file names')
    parser.add_argument('--annot-chr', dest='annot', required=True,
                        help='Annotation files for each chromosome')
    parser.add_argument('--frqfile-chr', dest='frqfile', required=True,
                        help='Frequency files for each chromosome')

    args = parser.parse_args(argv[1:])
    print(args)
    print("Loading gene annotations from '{}'".format(args.gene_ref))

    frq_file_df_list = []
    annot_file_df_list = []

    for chrom in range(1, 23):
        print(chrom)
        frq_file_partial_df = pd.read_csv("{}{}.frq".format(args.frqfile, chrom), delim_whitespace=True)
        maf_pass = (frq_file_partial_df.MAF > 0.05) & (frq_file_partial_df.MAF < 0.95)
        frq_file_df_list.append(frq_file_partial_df.loc[
                                    maf_pass])

        annot_file_partial_df = pd.read_csv("{}{}.annot.gz".format(args.annot, chrom), sep="\t")
        annot_file_df_list.append(annot_file_partial_df.loc[
                                      maf_pass])

    frq_file_df = pd.concat(frq_file_df_list, axis=0).reset_index(drop=True, inplace=False)
    annot_file_df = pd.concat(annot_file_df_list, axis=0).reset_index(drop=True, inplace=False)

    frq_file_df["BP"] = annot_file_df["BP"]

    annot_columns = annot_file_df.columns.drop(["CHR", "BP", "SNP", "CM"])

    for bed_file in args.bed_inclusion_file:
        inclusion_bed_file = pd.read_csv(bed_file, sep="\t", header=None, names=["chromosome", "start", "end", "name"])
        heritability_snps = count_heritability_snps_including(
            inclusion_bed_file, annot_columns, annot_file_df, frq_file_df)
        heritability_snps.T.map(str).agg(','.join, axis=1).to_csv(
            "M_5_50.{}.txt".format(bed_file.replace(".bed", "")), sep="\t")

    for bed_file in args.bed_exclusion_file:
        exclusion_bed_file = pd.read_csv(bed_file, sep="\t", header=None, names=["chromosome", "start", "end", "name"])
        heritability_snps = count_heritability_snps_excluding(
            exclusion_bed_file, annot_columns, annot_file_df, frq_file_df)
        heritability_snps.T.map(str).agg(','.join, axis=1).to_csv(
            "M_5_50.{}.txt".format(bed_file.replace(".bed", "")), sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main())
