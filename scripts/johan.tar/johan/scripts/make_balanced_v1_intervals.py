from gpn.data import load_table, get_balanced_intervals
import sys
_, parquet_file, gff_file, output_file, window_size, promoter_upstream = sys.argv

defined_intervals = load_table(parquet_file)
annotation = load_table(gff_file)
intervals = get_balanced_intervals(
    defined_intervals, annotation, window_size,
    promoter_upstream,
)
intervals.to_parquet(output_file, index=False)
