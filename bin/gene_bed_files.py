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

# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = argparse.ArgumentParser(description='Annotate significant eQTL results with variant and gene information')

    parser.add_argument('--gene-ref', dest='gene_ref', required=True,
                        help='Path to a GTF or GFF3 file')
    parser.add_argument('--out-prefix', dest='out_prefix', required=True,
                        help='Prefix to use for output file names')

    args = parser.parse_args(argv[1:])
    print(args)
    print("Loading gene annotations from '{}'".format(args.gene_ref))

    gencode_parser = GtfParser(args.gene_ref)
    gene_dataframe = gencode_parser.df

    assert np.alltrue(gene_dataframe.end > gene_dataframe.start)

    # In the dataframe, add a column to indicate what the right boundary (downstream) is for the trans window,
    # and what the left boundary (upstream) is for the trans window
    gene_dataframe["trans_downstream"] = gene_dataframe.end + 5*10**6
    gene_dataframe["trans_upstream"] = gene_dataframe.start - 5*10**6

    # In the dataframe, add a column to indicate what the right boundary (downstream) is for the cis window,
    # and what the left boundary (upstream) is for the cis window
    gene_dataframe["cis_downstream"] = gene_dataframe.end + 1*10**6
    gene_dataframe["cis_upstream"] = gene_dataframe.start - 1*10**6

    # Fix the upstream, lower boundary for both cis and the trans window to a minimum of 0
    gene_dataframe.loc[gene_dataframe["cis_upstream"] < 0, "cis_upstream"] = 0
    gene_dataframe.loc[gene_dataframe["trans_upstream"] < 0, "trans_upstream"] = 0

    print(gene_dataframe)

    # Output the boundaries.
    (gene_dataframe[["chromosome" , "trans_upstream", "trans_downstream", "gene_id"]]
     .to_csv("{}.trans.bed".format(args.out_prefix), sep="\t", index=False, header=False))
    (gene_dataframe[["chromosome" , "cis_upstream", "cis_downstream", "gene_id"]]
     .to_csv("{}.cis.bed".format(args.out_prefix), sep="\t", index=False, header=False))

    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
