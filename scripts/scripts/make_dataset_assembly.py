
from gpn.data import Genome, make_windows, get_seq
import math
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
import json

import sys
(_, parquet, fasta, *outputs, assembly, split_proportion, window_size, step_size, add_rc,
  whitelist_validation_chroms,whitelist_test_chroms) = sys.argv

splits = ['train', 'validation', 'test']


# Write split_proportion to a temporary file for debugging
debug_file_path = "/faststorage/project/johan_gpn/people/johanulsrup/johan_gpn/data/steps/tmp/split_proportion_debug.txt"
with open(debug_file_path, "w") as f:
    f.write(split_proportion)
print(f"split_proportion written to {debug_file_path}")

# Read the content of the temporary file to ensure it is correct
with open(debug_file_path, "r") as f:
    split_proportion_content = f.read()
    print(f"split_proportion_content: {split_proportion_content}")


split_proportion = json.loads(split_proportion)  # Parse the JSON string
split_proportions = [split_proportion[split] for split in splits]
assert np.isclose(sum(split_proportions), 1)

intervals = pd.read_parquet(input[0])
genome = Genome(input[1])
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
        assert os.path.isdir(dirname)
        os.makedir(dirname)
        print(dirname, 'created')
    print(path, split, dirname, os.getcwd(), os.path.exists(dirname))
    # to parquet to be able to load faster later
    print(intervals_split)
    print(intervals[(intervals_split==split).values])
    intervals[(intervals_split==split).values].to_parquet(
        path, index=False,
    )
