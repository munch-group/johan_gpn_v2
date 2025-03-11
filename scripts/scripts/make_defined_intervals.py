from gpn.data import Genome, filter_length
import sys

_, input_file, output_file, window_size = sys.argv

I = Genome(input_file).get_defined_intervals()
I = filter_length(I, window_size)
I.to_parquet(output_file, index=False)