
from gpn.data import (
    Genome, load_table, get_balanced_intervals, filter_length,
    filter_annotation_features,
)
import sys
import pandas as pd
_, parquet_file, fasta_file, output_file, window_size, feature, include_flank, add_jiter = sys.argv

I = pd.read_parquet(parquet_file)
annotation = load_table(fasta_file)
I = filter_annotation_features(
    I, annotation, feature,
    include_flank=include_flank, jitter=add_jitter,
)
I = filter_length(I, window_size)
I.to_parquet(output_file, index=False)
