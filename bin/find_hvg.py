import sys
import spatialdata as sp
import scanpy as sc
import pandas as pd
import argparse
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

def find_highly_variable_genes(norm_count_csv, output_plot, output_csv, is_scTransform):
    """
    This function takes as input a .csv file with normalized count data (NOT TRANSFORMED), creates an AnnData object, and finds highly variable genes 
    using the scanpy library. The specific parameters for the highly variable genes function are set to the Seurat defaults.

    This function saves a plot of the highly variable genes as a .png file.

    The plot to the left represent the non-centered dispersion values for each gene. 
    The plot to the right represents the mean centered, std scaled dispersion values for each gene.

    The values are reported in the log-scale

    This function also saves a .csv file with the list of genes, their mean, dispersions and dispersions_norm

    :param norm_count_csv: Path to the .csv file with the normalized counts data
    :param output_plot: Path to the output directory for the .png image
    :param output_csv: Path to the output directory for the csv file
    :param is_scTransform: Whether the data has been transformed with scTransform
    :return: None     
    """

    filename = norm_count_csv.split('/')[-1].split('.')[0]
    method = filename.split('_')[-1]

    # Read csv with the normalized counts data
    df = pd.read_csv(norm_count_csv, index_col=0)
    # Drop cells with nans or infs
    df = df.dropna(axis = 0, how = 'any')
    df = df.replace([np.inf, -np.inf], np.nan).dropna(axis = 0)
    
    if method == 'scTransform' or method == 'spanorm-logpac' or method == 'spanorm-pearson' or method == 'spanorm':
        print(f"Data normalized with {method}, skipping log1p transformation.")
        adata = sc.AnnData(df.to_numpy())
    else:
        print(f"Data normalized with {method}, applying log1p transformation.")
        adata = sc.AnnData(np.log1p(df.to_numpy()))

    adata.obs_names = df.index.values
    adata.var_names = df.columns.values

    # Highly Variable Genes
    if method != 'scTransform' and method != 'spanorm-pearson':
        adata = sc.pp.highly_variable_genes(adata, flavor = 'seurat', inplace = False) # Expect log1p transformed data, so for consistency we have to do it

        # Plot and save
        # sc.pl.highly_variable_genes(adata, show = False, log = True)
        # plt.savefig(output_plot + '/' + filename + '_hvg_logscale.png')
        # sc.pl.highly_variable_genes(adata, show = False, log = False)
        # plt.savefig(output_plot + '/' + filename + '_hvg_noscale.png')

        # Save dataframe containing the list of genes, their mean, dispersions and dispersions_norm
        # hvg_df = pd.DataFrame(adata.var)
        adata.rename(columns = {'means':'means', 'dispersions': 'dispersions'}, inplace = True)
        adata.to_csv(output_csv + '/' + filename + '_hvg.csv')

    else:      
        # Save a .csv file with the same format as hvg_df
        hvg_df = pd.DataFrame(adata.var)
        hvg_df['highly_variable'] = False
        hvg_df['means'] = adata.X.mean(axis = 0)
        hvg_df['dispersions'] = adata.X.var(axis = 0)
        hvg_df['dispersions_norm'] = None

        hvg_df.to_csv(output_csv + '/' + filename + '_hvg.csv')


    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find highly variable genes')
    parser.add_argument('--norm_counts_csv', type=str, help='Path to the .csv file with the normalized counts data')
    parser.add_argument('--output_plot', type=str, help='Path to the output directory for the .png image', default = '.')
    parser.add_argument('--output_csv', type=str, help='Path to the output directory for the csv file')
    parser.add_argument('--is_scTransform', help='Whether the data has been transformed with scTransform', action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    norm_counts_csv = args.norm_counts_csv
    output_plot = args.output_plot
    output_csv = args.output_csv
    is_scTransform = args.is_scTransform

    find_highly_variable_genes(norm_counts_csv, output_plot, output_csv, is_scTransform)