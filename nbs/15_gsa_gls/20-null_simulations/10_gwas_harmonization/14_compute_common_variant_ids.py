import pickle
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse

import numpy as np
import pandas as pd


def _get_gwas_variants(f):
    """
    read all GWAS and find a set of common panel_variant_id_values
    """
    gwas_data = pd.read_table(f, usecols=["panel_variant_id", "zscore"])
    assert gwas_data["panel_variant_id"].is_unique
    assert gwas_data.shape == gwas_data.dropna().shape
    return f.name, set(gwas_data["panel_variant_id"])

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-gwas-dir",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--input-gwas-file-pattern",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--n-samples",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--n-jobs",
        required=False,
        default=1,
        type=int,
    )
    parser.add_argument(
        "--random-seed",
        required=False,
        default=0,
        type=int,
    )

    args = parser.parse_args()
    
    input_gwas_dir = Path(args.input_gwas_dir).resolve()
    assert input_gwas_dir.exists(), input_gwas_dir
    
    input_files = sorted(list(input_gwas_dir.glob(args.input_gwas_file_pattern)))
    len(input_files)
    
    # sample files
    np.random.seed(args.random_seed)
    input_files = np.random.choice(input_files, size=args.n_samples, replace=False)
    len(input_files)
    
    common_variants = set()
    last_n_var_ids = -1
    with ProcessPoolExecutor(max_workers=args.n_jobs) as executor:
        for gwas_file_name, gwas_variants in executor.map(_get_gwas_variants, input_files, chunksize=10):
            if len(common_variants) == 0:
                common_variants = gwas_variants
            else:
                common_variants = common_variants.intersection(gwas_variants)
    
            n_var_ids = len(common_variants)
            same_previous = n_var_ids == last_n_var_ids
            last_n_var_ids = n_var_ids
            print(
                f"{gwas_file_name}, # common variants: {n_var_ids} (same? {same_previous})",
                flush=True
            )
    
    with open(input_gwas_dir / "common_variant_ids.pkl", 'wb') as f:
        pickle.dump(common_variants, f, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    run()
