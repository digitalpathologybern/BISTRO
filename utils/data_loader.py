"""
BISTRO utilities: data loading, layer management, and FOV assignment.

Consolidates shared I/O logic from NormAnalyzerPipeline, TransAnalyzerPipeline,
and 03_HVG_benchmark.py into standalone functions usable from any evaluation script.

Coordinate convention
---------------------
All spatial coordinates are in **micrometers (um)**. Column names:
  - x_local_um, y_local_um   : cell centroid within its FOV
  - x_global_um, y_global_um : cell centroid in slide coordinates
  - fov_center_x_um, fov_center_y_um : FOV center in slide coordinates
  - area_um2                  : cell area in square micrometers
"""

import gc
import os

import numpy as np
import pandas as pd
import spatialdata as sp
from glob import glob
from scipy.sparse import issparse

import geopandas as gpd
from shapely.geometry import Point
from shapely.strtree import STRtree
import ast

from helpers import get_fov_size, rasterize_fov_grid

from shapely.affinity import affine_transform, scale as shapely_scale


# ============================================================================
# SpatialData loading
# Source: NormAnalyzerPipeline.load_spatialdata_object (norm_analyzer.py, line 260)
# ============================================================================

def load_spatialdata(zarr_path):
    """
    Read a SpatialData zarr archive and return the object.

    Handles the EntityID index convention used by MERSCOPE datasets.

    Parameters
    ----------
    zarr_path : str
        Path to the .zarr directory.

    Returns
    -------
    spatialdata.SpatialData
        The loaded SpatialData object.
    """
    sd_obj = sp.SpatialData.read(zarr_path)

    if 'EntityID' in sd_obj.tables['filtered'].obs.columns:
        sd_obj.tables['filtered'].obs.set_index('EntityID', inplace=True, drop=False)
        sd_obj.tables['filtered'].obs.index.name = None
        if 'table' in sd_obj.tables:
            sd_obj.tables['table'].obs.set_index('EntityID', inplace=True, drop=False)
            sd_obj.tables['table'].obs.index.name = None

    print(f"SpatialData loaded: {sd_obj.tables['filtered'].shape[0]} cells, "
          f"{sd_obj.tables['filtered'].shape[1]} genes")

    return sd_obj


# ============================================================================
# Normalization layer discovery and I/O
# Source: discover_norm_layers from 03_HVG_benchmark.py (line 198)
#         load_single_layer from 03_HVG_benchmark.py (line 232)
#         unload_layer from NormAnalyzerPipeline (norm_analyzer.py, line 324)
# ============================================================================

def discover_norm_layers(nextflow_path):
    """
    Scan the Nextflow output ``norm/`` directory and return a mapping
    of {layer_name: filepath} for each normalization CSV.

    Filters out metadata files, size factor files, and raw count files.

    Parameters
    ----------
    nextflow_path : str
        Root path of the Nextflow output directory.

    Returns
    -------
    dict
        Mapping of normalization method name to its CSV file path.
    """
    norm_dir = os.path.join(nextflow_path, "norm")
    if not os.path.isdir(norm_dir):
        print(f"Warning: norm directory not found at {norm_dir}")
        return {}

    norm_map = {}
    for fpath in glob(os.path.join(norm_dir, "*.csv")):
        basename = os.path.basename(fpath)
        if any(tok in basename for tok in ("metadata", "SF", "counts.csv")):
            continue
        layer_name = os.path.splitext(basename)[0].split("_")[-1]
        norm_map[layer_name] = fpath

    print(f"Discovered {len(norm_map)} normalization layers: {sorted(norm_map.keys())}")
    return norm_map


def load_layer(adata, layer_name, norm_map):
    """
    Read a normalization CSV into adata.layers[layer_name] as float32.

    If the layer is already loaded or the name is 'X', this is a no-op.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to load the layer into.
    layer_name : str
        Name of the normalization method (key in norm_map).
    norm_map : dict
        Mapping of layer names to file paths.
    """
    if layer_name == "X" or layer_name in adata.layers:
        return

    if layer_name not in norm_map:
        raise ValueError(f"Layer '{layer_name}' not found in norm_map. "
                         f"Available: {sorted(norm_map.keys())}")

    print(f"Loading layer: {layer_name}")
    df = pd.read_csv(norm_map[layer_name], index_col=0).astype("float32")
    adata.layers[layer_name] = df.values
    del df
    gc.collect()


def load_layer_into_spatialdata(sd_obj, layer_name, norm_map):
    """
    Read a normalization CSV into sd_obj.tables['filtered'].layers[layer_name].

    Used by evaluation modules that operate on the SpatialData object directly
    (e.g., batch_effect.py which needs obs metadata alongside expression).

    Parameters
    ----------
    sd_obj : spatialdata.SpatialData
        The SpatialData object.
    layer_name : str
        Name of the normalization method.
    norm_map : dict
        Mapping of layer names to file paths.
    """
    if layer_name in sd_obj.tables['filtered'].layers:
        return

    if layer_name not in norm_map:
        raise ValueError(f"Layer '{layer_name}' not found in norm_map.")

    print(f"Loading layer into SpatialData: {layer_name}")
    df = pd.read_csv(norm_map[layer_name], index_col=0).astype("float32")
    sd_obj.tables['filtered'].layers[layer_name] = df
    del df
    gc.collect()


def unload_layer(obj, layer_name):
    """
    Remove a layer from memory and force garbage collection.

    Works with both AnnData and SpatialData objects.

    Parameters
    ----------
    obj : anndata.AnnData or spatialdata.SpatialData
        The object from which to remove the layer.
    layer_name : str
        Name of the layer to unload.
    """
    # Determine the target AnnData
    if hasattr(obj, 'tables'):
        adata = obj.tables['filtered']
    else:
        adata = obj

    if layer_name in adata.layers:
        print(f"Unloading layer: {layer_name}")
        del adata.layers[layer_name]

    # Also clean up any derived PCA or uns data for this layer
    keys_to_remove = [k for k in adata.uns.keys() if layer_name in k]
    for k in keys_to_remove:
        del adata.uns[k]

    keys_to_remove = [k for k in adata.varm.keys() if layer_name in k]
    for k in keys_to_remove:
        del adata.varm[k]

    if hasattr(adata, 'obsm'):
        keys_to_remove = [k for k in adata.obsm.keys() if layer_name in k]
        for k in keys_to_remove:
            del adata.obsm[k]

    gc.collect()


# ============================================================================
# Dense layer extraction
# Source: get_dense_layer from 03_HVG_benchmark.py (line 246)
# ============================================================================

def get_dense_layer(adata, layer_name):
    """
    Return an expression matrix for the specified layer as a dense float32 array.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object.
    layer_name : str
        Name of the layer, or 'X' for the default matrix.

    Returns
    -------
    np.ndarray
        Dense float32 expression matrix (cells x genes).
    """
    X = adata.X if layer_name == "X" else adata.layers[layer_name]
    if issparse(X):
        return np.asarray(X.toarray(), dtype=np.float32)
    return np.asarray(X, dtype=np.float32)


# ============================================================================
# FOV loading from enriched metadata
# ============================================================================

def load_fov_from_metadata(sd_obj, nextflow_path):
    """
    Load FOV assignments from the enriched metadata CSV (produced by
    bin/preprocessing/assign_fov.py) and merge them into the SpatialData
    obs table.

    This is the preferred method for obtaining FOV identifiers in
    evaluation scripts. The enriched metadata CSV contains 'fov',
    'fov_center_x_um', and 'fov_center_y_um' columns regardless of the
    iST platform.

    Parameters
    ----------
    sd_obj : spatialdata.SpatialData
        The SpatialData object whose obs table will be updated.
    nextflow_path : str
        Root path of the Nextflow output directory (containing norm/).

    Raises
    ------
    FileNotFoundError
        If no metadata CSV is found in the norm/ subdirectory.
    """
    norm_dir = os.path.join(nextflow_path, "norm")
    meta_files = glob(os.path.join(norm_dir, "*_metadata.csv"))
    if not meta_files:
        raise FileNotFoundError(
            f"No metadata CSV found in {norm_dir}. "
            f"Ensure readZarr and assignFOV have run."
        )

    meta = pd.read_csv(meta_files[0], index_col=0)
    adata = sd_obj.tables['filtered']

    fov_cols = [c for c in ['fov', 'fov_center_x_um', 'fov_center_y_um'] if c in meta.columns]
    if not fov_cols:
        print("Warning: no FOV columns found in metadata CSV.")
        return

    # Align metadata index with adata obs index
    meta.index = meta.index.astype(str)
    adata.obs.index = adata.obs.index.astype(str)

    for col in fov_cols:
        adata.obs[col] = meta.loc[adata.obs.index, col].values

    print(f"Loaded FOV from metadata: {fov_cols}, "
          f"{adata.obs['fov'].nunique()} unique FOVs")


def load_fov_info(nextflow_path):
    """
    Load the FOV info JSON produced by assign_fov.py.

    Parameters
    ----------
    nextflow_path : str
        Root path of the Nextflow output directory.

    Returns
    -------
    dict or None
        FOV info dictionary with keys: fov_source, technology, n_fovs, message.
        Returns None if the file is not found.
    """
    import json
    info_path = os.path.join(nextflow_path, "norm", "fov_info.json")
    if not os.path.exists(info_path):
        return None
    with open(info_path, 'r') as f:
        return json.load(f)


def assign_fov(sd_obj, technology, scale=1):
    """
    .. deprecated::
        Use ``bin/preprocessing/assign_fov.py`` as a Nextflow process instead.
        FOV assignment is now centralized and runs once after readZarr.
        Use ``load_fov_from_metadata()`` to read the pre-computed FOV
        identifiers in evaluation scripts.
    """
    import warnings
    warnings.warn(
        "assign_fov() is deprecated. FOV assignment is now centralized in "
        "bin/preprocessing/assign_fov.py. Use load_fov_from_metadata() instead.",
        DeprecationWarning, stacklevel=2
    )

    if technology == 'CosMx':
        print("Technology: CosMx (native FOV information available)")
        return

    fov_size = get_fov_size(technology, scale)
    print(f"Rasterizing pseudo-FOVs for {technology} with tile size {fov_size} um")

    sd_obj.tables['filtered'].obs = rasterize_fov_grid(
        sd_obj.tables['filtered'].obs, fov_size
    )
    if 'table' in sd_obj.tables:
        sd_obj.tables['table'].obs = rasterize_fov_grid(
            sd_obj.tables['table'].obs, fov_size
        )


# ============================================================================
# Tissue annotation loading
# ============================================================================

def load_tissue_annotations(sd_obj, annotation_path, he_alignment_path=None,
                            pixel_size_um=0.2125):
    """
    Load tissue region annotations from a CSV or GeoJSON file and merge
    them into the SpatialData obs table.

    For CSV files: expects a column 'tissue_annotations' (or 'niche' /
    'banksy_0.8' which are renamed) and either 'cell_ID' or 'Unnamed: 0'
    as the merge key.

    For GeoJSON files: optionally applies an H&E-to-Xenium alignment affine,
    then scales polygons from Xenium pixels to um, then performs
    point-in-polygon assignment against x_global_um / y_global_um via
    Shapely STRtree.

    Parameters
    ----------
    sd_obj : spatialdata.SpatialData
    annotation_path : str
        Path to the annotation file (.csv or .geojson).
    he_alignment_path : str, optional
        Path to a 3x3 CSV (no header) containing the H&E-to-Xenium affine.
        Applied to polygons before scaling to um. Only used for GeoJSON.
    pixel_size_um : float, optional
        Pixel size in um, used to scale post-affine polygons to um.
        Default 0.2125.
    """
    if not annotation_path or not os.path.exists(annotation_path):
        print(f"Warning: annotation path not found: {annotation_path}")
        return

    ext = annotation_path.rsplit('.', 1)[-1].lower()

    if ext == 'csv':
        anno_df = pd.read_csv(annotation_path)

        # Normalize column names
        if "niche" in anno_df.columns:
            anno_df.rename(columns={"niche": "tissue_annotations"}, inplace=True)
        if "banksy_0.8" in anno_df.columns:
            anno_df.rename(columns={"banksy_0.8": "tissue_annotations"}, inplace=True)
            anno_df = anno_df[['Unnamed: 0', 'tissue_annotations']]
            anno_df['tissue_annotations'] = anno_df['tissue_annotations'].astype(str)

        # Build merge key
        if 'cell_ID' in anno_df.columns:
            obs_index = set(sd_obj.tables['filtered'].obs.index.astype(str))
            if anno_df['cell_ID'].astype(str).iloc[0] in obs_index:
                anno_df['merge_key'] = anno_df['cell_ID'].astype(str)
            else:
                anno_df['merge_key'] = anno_df['cell_ID'].apply(
                    lambda x: f"{x.split('_')[-1]}_{x.split('_')[-2]}"
                )
        elif 'Unnamed: 0' in anno_df.columns:
            anno_df['merge_key'] = anno_df['Unnamed: 0']
        else:
            print("Warning: cannot determine merge key in annotation CSV.")
            return

        # Filter to matching cells
        obs_index = list(sd_obj.tables['filtered'].obs.index.astype(str))
        anno_df = anno_df[anno_df['merge_key'].astype(str).isin(obs_index)]
        anno_df = anno_df.set_index('merge_key')
        if 'Unnamed: 0' in anno_df.columns:
            anno_df.index.name = None

        # Merge into obs
        sd_obj.tables['filtered'].obs.index = sd_obj.tables['filtered'].obs.index.astype(str)
        anno_df.index = anno_df.index.astype(str)

        drop_cols = [c for c in ['cell_ID', 'Unnamed: 0'] if c in anno_df.columns]
        sd_obj.tables['filtered'].obs = sd_obj.tables['filtered'].obs.join(
            anno_df.drop(columns=drop_cols, errors='ignore'), how='left'
        )

        n_annotated = sd_obj.tables['filtered'].obs['tissue_annotations'].notna().sum()
        print(f"Tissue annotations loaded (CSV): {n_annotated} cells annotated")

    elif ext == 'geojson':
        def _extract_name(classification):
            if isinstance(classification, dict):
                return classification.get("name")
            if isinstance(classification, str):
                try:
                    parsed = ast.literal_eval(classification)
                    if isinstance(parsed, dict):
                        return parsed.get("name")
                except (ValueError, SyntaxError):
                    return None
            return None

        # Load affine if provided
        affine_params = None
        if he_alignment_path and os.path.exists(he_alignment_path):
            M = pd.read_csv(he_alignment_path, header=None).values
            if M.shape != (3, 3):
                raise ValueError(f"Expected 3x3 affine, got {M.shape} from {he_alignment_path}")
            affine_params = [M[0][0], M[0][1], M[1][0], M[1][1], M[0][2], M[1][2]]
            print(f"Loaded alignment affine from {he_alignment_path}")

        anno_df = gpd.read_file(annotation_path)
        anno_df['category'] = anno_df["classification"].apply(_extract_name)

        # Apply affine + scale to um. 10x Xenium alignment CSVs map H&E px
        # to Xenium px, so after the affine we scale by pixel_size.
        if affine_params is not None:
            anno_df["transformed_geometry"] = anno_df["geometry"].apply(
                lambda g: shapely_scale(
                    affine_transform(g, affine_params),
                    xfact=pixel_size_um,
                    yfact=pixel_size_um,
                    origin=(0, 0),
                )
            )
        else:
            anno_df["transformed_geometry"] = anno_df["geometry"]

        # Build spatial index
        geoms = list(anno_df["transformed_geometry"])
        cats = list(anno_df["category"])
        tree = STRtree(geoms)
        index_to_category = {i: cats[i] for i in range(len(geoms))}

        # Point-in-polygon assignment using global um coordinates
        x_coords = sd_obj.tables["filtered"].obs["x_global_um"].values
        y_coords = sd_obj.tables["filtered"].obs["y_global_um"].values

        cell_labels = []
        for x, y in zip(x_coords, y_coords):
            pt = Point(x, y)
            match_indices = tree.query(pt)
            found_label = None
            for idx in match_indices:
                if geoms[idx].contains(pt):
                    found_label = index_to_category[idx]
                    break
            cell_labels.append(found_label if found_label else "None")

        sd_obj.tables["filtered"].obs["tissue_annotations"] = cell_labels
        n_annotated = sum(1 for l in cell_labels if l != "None")
        print(f"Tissue annotations loaded (GeoJSON): {n_annotated} cells annotated")

    else:
        print(f"Warning: unrecognized annotation format '{ext}'. Use .csv or .geojson.")

# ============================================================================
# Spatial coordinate extraction
# Source: get_spatial_coordinates from 03_HVG_benchmark.py (line 254)
# ============================================================================

def get_spatial_coordinates(adata, sd_obj=None):
    """
    Extract 2-D spatial coordinates aligned with adata.obs_names.

    Tries multiple common column name conventions across iST platforms,
    prioritizing the standardized um-unit columns.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with cell metadata in obs.
    sd_obj : spatialdata.SpatialData, optional
        If provided, also checks shape centroids as a fallback.

    Returns
    -------
    np.ndarray
        Array of shape (n_cells, 2) with spatial coordinates in um.

    Raises
    ------
    ValueError
        If no spatial coordinates can be extracted.
    """
    if "spatial" in adata.obsm:
        return np.asarray(adata.obsm["spatial"])[:, :2].astype(np.float64)

    col_pairs = [
        ("x_global_um", "y_global_um"),
        ("x_local_um", "y_local_um"),
        ("x_centroid", "y_centroid"),
        ("CenterX_global_px", "CenterY_global_px"),   # legacy fallback
        ("CenterX_local_px", "CenterY_local_px"),      # legacy fallback
        ("x_location", "y_location"),
        ("center_x", "center_y"),
    ]
    for xc, yc in col_pairs:
        if xc in adata.obs.columns and yc in adata.obs.columns:
            return adata.obs[[xc, yc]].values.astype(np.float64)

    if sd_obj is not None:
        for shape_key in sd_obj.shapes:
            shapes = sd_obj.shapes[shape_key]
            if hasattr(shapes, "geometry"):
                centroids = shapes.geometry.centroid
                coords = np.column_stack([centroids.x, centroids.y])
                if coords.shape[0] == adata.n_obs:
                    return coords.astype(np.float64)

    raise ValueError("Could not extract spatial coordinates from any source.")