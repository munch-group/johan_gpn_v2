from gpn.data import Genome, make_windows, get_seq
import math
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
import ast

import sys
(_, parquet, fasta, *outputs, assembly, split_proportion, window_size, step_size, add_rc,
  whitelist_validation_chroms, whitelist_test_chroms) = sys.argv

splits = ['train', 'validation', 'test']

window_size = int(window_size)
step_size = int(step_size)
add_rc = bool(add_rc)

# Hardcoded split_proportion dictionary
split_proportion = {'train': 0.7, 'validation': 0.2, 'test': 0.1}

split_proportions = [split_proportion[split] for split in splits]
assert np.isclose(sum(split_proportions), 1)

intervals = pd.read_parquet(parquet)
genome = Genome(fasta)
intervals = make_windows(
    intervals, window_size, step_size, add_rc,
)

intervals = intervals.sample(frac=1.0, random_state=42)
intervals["assembly"] = assembly
intervals = intervals[["assembly", "chrom", "start", "end", "strand"]]
intervals = get_seq(intervals, genome)

chroms = intervals.chrom.unique()
chrom_split = np.random.choice(
    splits, p=split_proportions, size=len(chroms),
)
chrom_split[np.isin(chroms, whitelist_validation_chroms)] = "validation"
chrom_split[np.isin(chroms, whitelist_test_chroms)] = "test"
chrom_split = pd.Series(chrom_split, index=chroms)

intervals_split = chrom_split[intervals.chrom]

for path, split in zip(outputs, splits):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        print(dirname, 'created')
    print(path, split, dirname, os.getcwd(), os.path.exists(dirname))
    # to parquet to be able to load faster later
    print(intervals_split)
    print(intervals[(intervals_split==split).values])
    intervals[(intervals_split==split).values].to_parquet(
        path, index=False,
    )