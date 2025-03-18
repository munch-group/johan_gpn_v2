import pandas as pd
from tqdm import tqdm
import os
import math
from pathlib import Path
import sys

inputs = sys.argv[1:-1]  # Input paths are passed as arguments
output_dir = sys.argv[-1]  # Output directory is the last argument

config = {
    "subsample_to_target": True,
    "target_assembly": "some_assembly",
    "samples_per_file": 1000
}

try:
    intervals = pd.concat(
        tqdm([pd.read_parquet(path) for path in inputs], total=len(inputs)),
        ignore_index=True,
    ).sample(frac=1, random_state=42)
    print(intervals)
except Exception as e:
    print(f"Error reading parquet files: {e}")
    raise

if config.get("subsample_to_target", False) and "train" in inputs:
    n_target = (intervals.assembly == config["target_assembly"]).sum()
    intervals = intervals.groupby("assembly").sample(
        n=n_target, random_state=42
    ).sample(frac=1, random_state=42)
    print("train", intervals.assembly.value_counts())
    print(intervals)

n_shards = math.ceil(len(intervals) / config["samples_per_file"])
if n_shards >= 10000:
    print(f"Warning: Number of shards ({n_shards}) exceeds the limit of 10000. Reducing the number of samples per file.")
    n_shards = 9999

os.makedirs(output_dir, exist_ok=True)
for i in tqdm(range(n_shards)):
    path = Path(output_dir) / f"shard_{i:05}.jsonl.zst"
    intervals.iloc[i::n_shards].to_json(
        path, orient="records", lines=True,
        compression={'method': 'zstd', 'threads': -1}
    )
