

import sys

_, input_file, output_file, window_size = sys.argv

from gpn.data import  Genome, filter_length

print(f"Running make_all_intervals with input_file={input_file}, output_file={output_file}, window_size={window_size}")

try:
    I = Genome(input_file).get_all_intervals()
    I = filter_length(I, int(window_size))
    I.to_parquet(output_file, index=False)
    print(f"Intervals saved to {output_file}")
except Exception as e:
    print(f"Error: {e}")