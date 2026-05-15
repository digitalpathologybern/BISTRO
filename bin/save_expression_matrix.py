import sys
import spatialdata as sp
from scipy.sparse import issparse
import pandas as pd

def read_zarr_to_csv(zarr_path, output_csv_path, store_table):
    """ Reads a .zarr file and generate a .csv file of the gene expression matrix specified under "store_table".

    The generated .csv file has the shape cells by genes, where rows are numerated and column names are genes

    :param zarr_path: Path to zarr file
    :param output_csv_path: Path to output csv file. Do not include the name of the file
    :param store_table: Name of the table to be stored
    :return: None     
    """
    # Load the .zarr file
    z = sp.SpatialData.read(zarr_path)
    if 'EntityID' in z.tables['filtered'].obs.columns:
        z.tables['filtered'].obs.set_index('EntityID', inplace=True, drop = False)
        z.tables['table'].obs.set_index('EntityID', inplace=True, drop = False)
        z.tables['filtered'].obs.index.name = None
        z.tables['table'].obs.index.name = None
        
    # Get the name of the zarr file
    zarr_name = zarr_path.split('/')[-1][:-5]
    # Assuming you have the count matrix stored in 'X' of the AnnData object
    X = z.tables[store_table].X
    if issparse(X):
        count_matrix = pd.DataFrame(X.toarray())
    else:
        count_matrix = pd.DataFrame(X)
    count_matrix.columns = z.tables[store_table].var_names
    count_matrix.index = z.tables[store_table].obs.index

    # Get the metadata, where the area values are stored
    metadata_matrix = z.tables[store_table].obs

    # Save the count matrix to CSV
    count_matrix.to_csv(output_csv_path + '/' + zarr_name + '_' + store_table + '_counts.csv')
    metadata_matrix.to_csv(output_csv_path + '/' + zarr_name + '_' + store_table + '_metadata.csv')

if __name__ == "__main__":
    zarr_path = sys.argv[1]
    output_csv_path = sys.argv[2]
    store_table = sys.argv[3]
    read_zarr_to_csv(zarr_path, output_csv_path, store_table)