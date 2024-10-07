#!/usr/bin/env python3

import argparse
import os
import re
import sys
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import stats

SCHEMA = pa.schema([("variant_index", pa.int64()), ("beta", pa.float64()),
                    ("standard_error", pa.float64()), ("i_squared", pa.float64()),
                    ("sample_size", pa.float64())])


class QtlFilter(ABC):
    _operator = None
    _field_name = None
    @abstractmethod
    def _get_value(self):
        raise NotImplementedError
    def get_filter(self):
        return self._field_name, self._operator, self._get_value()
    def get_field_name(self):
        return self._field_name


class QtlSetFilter(QtlFilter):
    _operator = "in"
    def __init__(self, values):
        self._values = set(values)
    def _get_value(self):
        return self._values
    def get_values(self):
        return self._get_value()
    @classmethod
    def from_path(cls, path):
        df = pd.read_csv(path, header = None, delimiter = '\t')
        return cls(df.iloc[:,0].tolist())
    @classmethod
    def from_list(cls, genes):
        return cls(genes)
    def apply_to_df(self, df):
        return df.loc[df[self._field_name].isin(self._values)]


class QtlVariantFilter(QtlSetFilter):
    _field_name = "variant_index"
    def __init__(self, variants):
        super().__init__(variants)


class QtlGeneFilter(QtlSetFilter):
    _field_name = "phenotype"
    def __init__(self, genes):
        super().__init__(genes)


class QtlCohortFilter(QtlSetFilter):
    _field_name = "cohort"
    def __init__(self, cohorts):
        super().__init__(cohorts)


class QtlChromosomeFilter(QtlSetFilter):
    _field_name = "chromosome"
    _operator = "=="
    def __init__(self, chromosome):
        self._value = chromosome
    def _get_value(self):
        return self._value


class QtlLocusVariantFilter:
    def __init__(self, variants):
        self._variant_filter = QtlVariantFilter(variants)
    def get_filters(self):
        # Return a list, with in each list a chromosome filter, and a
        return [self._variant_filter.get_filter()]
    @classmethod
    def from_locus(cls, chromosome, start, stop, variant_reference):
        chromosome_reference = variant_reference.loc[variant_reference["chromosome"] == chromosome, :]
        variants = chromosome_reference.loc[
            ((chromosome_reference["bp"] >= start) & (chromosome_reference["bp"] <= stop)), "variant_index"]
        return cls(variants)


class QtlPThresholdFilter(QtlFilter):
    _field_name = "p_value"
    _operator = "<"
    def __init__(self, p_value):
        self._value = p_value
    def apply_to_df(self, df):
        return df.loc[df[self._field_name] < self._value]
    def apply(self, series):
        return series < self._value
    def _get_value(self):
        return self._value


class QtlNThresholdFilter(QtlFilter):
    _field_name = "sample_size"
    _operator = ">="
    def __init__(self, n_threshold):
        self._n = n_threshold
    def apply_to_df(self, df):
        return df.loc[df[self._field_name] >= self._n]
    def apply(self, series):
        return series >= self._n
    def _get_value(self):
        return self._n


class QtlResultProcessor:
    def __init__(self, path, gene_filter=None):
        self.path = path
        self.gene_filter = gene_filter
        self.variant_filters = None
        self.significance_filter = None
        self.n_filter = None
        self.cohort_filter = None
        self.column_mapping = \
            {"t_stat": self.t_stat,
             "z_score": self.z_score,
             "p_value": self.p_value
             }
        self._t_stat = None
        self._z_score = None
        self._p_value = None
        self._df = None
    def extract(self,
                cols=None,
                drop=None,
                add=None):
        if add is None:
            add = set()
        if drop is None:
            drop = set()
        if cols is None:
            cols = set()
        filters = self.get_filters()
        print(filters)
        dataset = pq.ParquetDataset(
            self.path, validate_schema=True,
            filters=filters)
        self._df = dataset.read().to_pandas()
        print("Queried {} rows".format(self._df.shape[0]))
        cols_to_add = (
            cols.union(add)
            .intersection(self.column_mapping.keys()))
        print("Columns to add: {}".format(cols_to_add))
        self.include_cols(cols_to_add)
        print("Added columns. Table now has {} columns".format(self._df.shape[1]))        
        if self.n_filter is not None:
            self._df = self._df.loc[self.n_filter.apply(self._df.sample_size)]
            print("Filtered on samle size. {} rows remain".format(self._df.shape[0]))
        if self.significance_filter is not None:
            self._df = self._df.loc[self.significance_filter.apply(self.p_value())]
            self._p_value = self._p_value[self.significance_filter.apply(self.p_value())]
            print("Filtered to {} rows".format(self._df.shape[0]))
        if cols is not None and len(cols) > 0:
            default_cols = ["variant_index", "phenotype"]
            default_cols.extend(cols.union(add))
            out = self._df.loc[:, default_cols]
        elif drop is not None and len(drop) > 0:
            out = self._df.drop(drop, axis=1)
        else:
            out = self._df
        return out
    def t_stat(self):
        if self._t_stat is None:
            self._t_stat = self.calculate_t_stat()
        return self._t_stat
    def p_value(self):
        if self._p_value is None:
            self._p_value = self.calculate_p_value()
        return self._p_value
    def z_score(self):
        if self._z_score is None:
            self._z_score = self.calculate_t_stat()
        return self._z_score
    def calculate_t_stat(self):
        return self._df['beta']/self._df['standard_error']
    def calculate_z_score(self):
        return self.calculate_t_stat()
    def calculate_p_value(self, degrees_of_freedom = None):
        if degrees_of_freedom is None:
            return stats.norm.sf(np.abs(self.z_score()))*2
        else:
            return stats.t.sf(np.abs(self.z_score()), degrees_of_freedom)*2
    def get_filters(self):
        base_filter_list = [self.gene_filter.get_filter()]
        if self.cohort_filter is not None:
            base_filter_list += [self.cohort_filter.get_filter()]
        print(base_filter_list)
        if self.variant_filters is not None:
            print(self.variant_filters)
            filter_list = [base_filter_list + [filter.get_filter()] for filter in self.variant_filters]
        else:
            filter_list = [base_filter_list]
        return filter_list
    def include_cols(self, cols):
        for col in sorted(cols):
            self._df[col] = self.column_mapping[col]()
    def add_n_threshold_filter(self, filter):
        self.n_filter = filter


def column_specification(cols):
    column_specifications = {"-": set(), "+": set(), "": set()}
    if cols is not None and cols != "":
        for column_specification in cols.split(","):
            regex_match = re.match(r"([+-]?)(\w+)", column_specification)
            if not regex_match:
                argparse.ArgumentTypeError(
                    "Column specification {} did not match expected pattern"
                    .format(column_specification))
            column_specifications[regex_match.group(1)].add(regex_match.group(2))
    return column_specifications


def export_write(input_file, output_prefix, qtl_gene_filter, variant_filter, cohort_filter, column_specifications, p_thresh=None, as_matrix=False):
    file_conns = dict()
    columns_to_write = column_specifications[""].union(column_specifications["+"])
    try:
        if as_matrix:
            for col in columns_to_write:
                file_conns[col] = open("{}.out.{}.csv".format(output_prefix, col), 'w')
            result_processor = QtlResultProcessor(
                input_file, qtl_gene_filter)
            print(variant_filter)
            result_processor.variant_filters = [variant_filter] if variant_filter is not None else None
            result_processor.cohort_filter = cohort_filter
            if p_thresh is not None:
                result_processor.significance_filter = QtlPThresholdFilter(p_thresh)
            df = result_processor.extract(
                cols=column_specifications[""],
                drop=column_specifications["-"],
                add=column_specifications["+"])
            for col in columns_to_write:
                pd.pivot(df, columns="phenotype", index="variant_index", values=col).to_csv(file_conns[col], sep="\t", header=True, index=True, na_rep="NA")
        else:
            first = True
            file_conns["long"] = open("{}.out.csv".format(output_prefix), 'w')
            for gene in qtl_gene_filter.get_values():
                print("Gene {}".format(gene))
                qtl_single_gene_filter = QtlGeneFilter.from_list([gene])
                result_processor = QtlResultProcessor(
                    input_file, qtl_single_gene_filter)
                result_processor.variant_filters = [variant_filter] if variant_filter is not None else None
                result_processor.cohort_filter = cohort_filter 
                if p_thresh is not None:
                    result_processor.significance_filter = QtlPThresholdFilter(p_thresh)
                df = result_processor.extract(
                    cols=column_specifications[""],
                    drop=column_specifications["-"],
                    add=column_specifications["+"])
                df.to_csv(file_conns["long"], sep="\t", header=first, index=None)
                first = False
    finally:
        print("Done!")
        for file_conn in file_conns.values():
            file_conn.close()
            print("Done!")
            print("Closing output file '{}'".format(os.path.basename(file_conn.name)))


def export_write_qtl_pairs(input_file, output_file, qtl_gene_variant_df, cohort_filter, column_specifications, p_thresh=None, as_matrix=False):
    first = True
    with open(output_file, 'w') as f:
        for gene,chunk in qtl_gene_variant_df.groupby('gene'):
            print("Gene {}".format(gene))
            print(chunk)
            qtl_single_gene_filter = QtlGeneFilter.from_list([gene])
            result_processor = QtlResultProcessor(
                input_file, qtl_single_gene_filter)
            variant_filters = [QtlVariantFilter(chunk['variant_index'])]
            result_processor.variant_filters = variant_filters
            result_processor.cohort_filter = cohort_filter
            if p_thresh is not None:
                result_processor.significance_filter = QtlPThresholdFilter(p_thresh)
            df = result_processor.extract(
                cols=column_specifications[""],
                drop=column_specifications["-"],
                add=column_specifications["+"])
            df.to_csv(f, sep="\t", header=first, index=None)
            first = False
        print("Done!")
        print("Closing output file '{}'".format(output_file))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description = "Analyze HASE output parquet files.")

    parser.add_argument('-i', '--input-file', required = True,
                        help = """One or multiple input parquet arrays. 
                                  Usage of wildcards is supported but then the argument has to be quoted.""")
    parser.add_argument('-o', '--output-prefix', type = str,
                        required = True,
                        help = "Path to tab-separated output file.")
    parser.add_argument('-p', '--p-thresh', type = float,
                        required = False,
                        help = "P-value threshold for filtering the results. Defaults None.",
                        default = None)
    parser.add_argument('-g', '--genes', required = False, default = None, nargs = '+',
                        help = """Individual phenotype IDs specified and separated by space.""")
    parser.add_argument('-G', '--genes-file', required = False, default = None,
                        help = """File with the list of phenotypes to include.""")
    parser.add_argument('-v', '--variants', required = False, default = None, nargs = '+',
                        help = """Individual SNP IDs specified and separated by space.""")
    parser.add_argument('-V', '--variants-file', required = False, default = None,
                        help = """File with the list of SNPs/variants to include.""")
    parser.add_argument('-B', '--bed-file', required=False, default=None,
                        help = """Bed file with a list of loci to extract""")
    parser.add_argument('-P', '--gene-variant-pair-file', required=False, default=None, 
                        help="""Gene variant pairs to extract""")
    parser.add_argument('-r', '--variant-reference', required = False,
                        help = "Reference for variants. Has to be gzipped and space-delimited.")
    parser.add_argument('-c', '--cols', dest="column_specifications", required = False, default = None,
                        type=column_specification, help="""Extract only z-scores""")
    parser.add_argument('-n', '--n-threshold', required=False, default=None,
                        help = "Minimal sample size")
    parser.add_argument('-m', '--as-matrix', dest="as_matrix", required=False, default=False, action='store_true')
    parser.add_argument('-C', '--cohort', required=False, default=None)

    args = parser.parse_args(argv[1:])

    print(args)

    variant_reference = None
    variants_list = None
    variant_filter = None
    cohort_filter = None
    loci = None

    if args.cohort is not None:
        cohort_filter = QtlCohortFilter({args.cohort,})

    if args.n_threshold is not None:
        n_threshold_filter = QtlNThresholdFilter(args.n_threshold)

    if args.variant_reference is not None and (args.variants is not None or args.variants_file is not None or args.bed_file is not None or args.gene_variant_pair_file is not None):
        variant_reference = pd.read_parquet(args.variant_reference)
        print(variant_reference.head())

    if args.gene_variant_pair_file is not None:
        gene_variant_pairs = pd.read_csv(args.gene_variant_pair_file, delimiter="\t", header=None, names=["variant", "gene"])
        variant_list = gene_variant_pairs["variant"]
        if variant_reference is None:
            parser.error("Cannot subset on variants without variant reference")
        processed_gene_variant_pairs = pd.merge(gene_variant_pairs, variant_reference, left_on="variant", right_on="variant")
        print("Starting export")
        output_file = "{}.out.csv".format(args.output_prefix)
        export_write_qtl_pairs(args.input_file, output_file,
                     processed_gene_variant_pairs, cohort_filter,
                     args.column_specifications, args.p_thresh)
        return 0

    if args.variants_file is not None:
        print("Using variants file '%s' to filter on variants." % args.variants_file)
        if variants_list is not None:
            print("Variant filter already defined. Skipping...")
        variants_list = (
            pd.read_csv(args.variants_file, header=None, delimiter='\t')
            .iloc[:,0].tolist())
    if args.variants is not None:
        print("Provided %d variants for filtering." % len(args.variants))
        if variants_list is not None:
            print("Variant filter already defined. Skipping...")
        variants_list = args.variants
    if args.bed_file is not None:
        print("Using loci file '%s' to filter on variants." % args.variants_file)
        if variants_list is not None:
            print("Variant filter already defined. Skipping...")
        loci = pd.read_csv(args.bed_file, sep="\t", header=None, names=["chromosome", "start", "stop", "name"])

    if variants_list is not None:
        if variant_reference is None:
            parser.error("Cannot subset on variants without variant reference")
        variant_selection = variant_reference.loc[variant_reference.loc[:, "variant"].isin(variants_list), "variant_index"]
        variant_filter = QtlVariantFilter(variant_selection)
        print("Will filter on %d variants" % len(variant_selection))

    qtl_gene_filter = None

    if args.genes_file is not None:
        print("Using genes file '%s' to filter on variants." % args.genes_file)
        qtl_gene_filter = QtlGeneFilter.from_path(args.genes_file)
    if args.genes is not None:
        print("Provided %d genes for filtering." % len(args.genes))
        if qtl_gene_filter is not None:
            print("Variant filter already defined. Skipping...")
        qtl_gene_filter = QtlGeneFilter.from_list(args.genes)

    if loci is None:
        print("Starting export")
        export_write(args.input_file, args.output_prefix,
                     qtl_gene_filter, variant_filter, cohort_filter,
                     args.column_specifications, args.p_thresh, args.as_matrix)

    else:
        for i, (index, row) in enumerate(loci.iterrows()):
            print("Starting export for locus {}/{}".format(i+1, loci.shape[0]))

            locus = row["name"].split(",")
            chromosome = row["chromosome"]
            start = row["start"]
            stop = row["stop"]

            locus_filter = QtlLocusVariantFilter.from_locus(
                chromosome, start, stop, variant_reference)

            if qtl_gene_filter is None:
                qtl_gene_filter = QtlGeneFilter.from_list(locus)

            output_file = "{}.{}_{}-{}.out.csv".format(args.output_prefix, chromosome, start, stop)

            export_write(args.input_file, output_file,
                         qtl_gene_filter, [locus_filter], cohort_filter,
                         args.column_specifications, args.p_thresh, args.as_matrix)

    return 0


if __name__ == '__main__':
    sys.exit(main())
