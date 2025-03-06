

import sys

_, input_file, output_file, window_size = sys.argv

from gpn.data import  Genome, filter_length

I = Genome(input_file).get_all_intervals()
I = filter_length(I, int(window_size))
I.to_parquet(output_file, index=False)
