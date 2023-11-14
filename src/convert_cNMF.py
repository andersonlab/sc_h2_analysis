#!/usr/bin/env python

__author__ = 'Tobi Alegbe'
__date__ = '2023-09-21'
__version__ = '0.0.2'


import argparse
import pandas as pd
import numpy as np

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read cNMF gene spectra file and convert to CELLECT input.
            """
    )

    parser.add_argument(
        '-cnmf', '--cnmf_spectra',
        action='store',
        dest='cnmf_spectra',
        required=True,
        help='Gene spectra output from cNMF'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='output_file',
        default='sc_CELLEX_data',
        help='Basename of output files.\
            (default: %(default)s)'
    )

    options = parser.parse_args()
    gene_score_file = options.cnmf_spectra
    output_file = options.output_file
    
    gene_scores = pd.read_csv(gene_score_file, sep='\t', index_col=0).T


    # Add text prefix to column names
    gene_scores = gene_scores.add_prefix('Factor_')
    gene_scores.index.name = 'ENSG'
    gene_scores

    # Normalise gene scores to min-0.99 quantile (not max as max is very large)
    norm_gene_scores = gene_scores/gene_scores.quantile(0.99, axis=0)
    norm_gene_scores = norm_gene_scores.clip(0,1)

    # Add text prefix to column names
    norm_gene_scores.index.name = 'gene'
    norm_gene_scores

    norm_gene_scores.to_csv(
                '{}.nmf.csv.gz'.format(output_file),
                compression='gzip',
                index=True,
                header=True
        )

if __name__ == '__main__':
    main()