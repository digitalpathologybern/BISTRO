import pandas as pd
import sys
import os

def normalize_with_area(raw_counts, raw_metadata, area_norm):
    """ Reads a .csv file containing raw RNA detections and save a .csv with the area-normalized values.

    The generated .csv file has the shape cells by genes, where rows are numerated and column names are genes

    :param raw_counts: Path to raw .csv file
    :param raw_metadata: Path to the metadata file of the raw counts
    :param area_norm: Path to output .csv file. Do not include the name of the file
    :return: None     
    """
    # Load the .csv files
    # print("Loading counts")
    # counts = pd.read_csv(raw_counts, index_col = 0)
    # print("counts loaded")
    print("Loading metadata")
    meta = pd.read_csv(raw_metadata, index_col = 0)
    print("metadata loaded")
    # Get the name of the zarr file
    counts_name = raw_counts.split('/')[-1][:-4]
    # Get the area of each cell
    area = meta['area_um2']
    median = area.median()
    area = area/median

    ###### This is memory intensive ######
    # norm = counts.div(area, axis = 'rows')
    # # Save the count matrix to CSV
    # norm.to_csv(area_norm + '/' + counts_name + '_areaNorm.csv')
    # area.to_csv(area_norm + '/' + counts_name + '_areaSF.csv')


    ##### Do this instead
    out_counts = os.path.join(area_norm, f"{counts_name}_areaNorm.csv")
    out_sf = os.path.join(area_norm, f"{counts_name}_areaSF.csv")

    # --- Write header with correct columns ---
    reader = pd.read_csv(raw_counts, index_col=0, chunksize=5000)
    first_chunk = next(reader)
    first_chunk = first_chunk.div(area[first_chunk.index], axis=0)
    first_chunk.to_csv(out_counts)
    del first_chunk
    # Continue with remaining chunks
    for chunk in reader:
        norm_chunk = chunk.div(area[chunk.index], axis=0)
        norm_chunk.to_csv(out_counts, mode='a', header=False)
        del norm_chunk

    # --- Save area size factors ---
    area.to_csv(out_sf)


if __name__ == "__main__":
    raw_counts = sys.argv[1]
    raw_metadata = sys.argv[2]
    area_norm = sys.argv[3]
    normalize_with_area(raw_counts, raw_metadata, area_norm)