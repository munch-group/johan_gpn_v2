from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import bioframe as bf
from collections import defaultdict
from gpn.data import (
    filter_defined, filter_length, load_table, add_flank, get_annotation_features,
    add_jitter, get_promoters, get_random_intervals, union_intervals,
    intersect_intervals, intervals_size
)
import gzip
import h5py
import math
#import more_itertools
import numpy as np
import os
import pandas as pd
import re
import scipy.sparse as sp_sparse
from scipy.special import softmax
from scipy.stats import entropy
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import torch
from tqdm import tqdm
tqdm.pandas()

from gpn.data import make_windows, get_seq
from gpn.data import load_fasta, save_fasta, load_table, load_repeatmasker, Genome

from gwf import Workflow
import yaml

gwf = Workflow()

configfile = "scripts/single_genom_analysis/config.yaml"  # this is the path to the config file

print("hej")

# # Load the configuration file
# with open(configfile, 'r') as file:
#     config = yaml.safe_load(file)

# # Ensure config is defined before accessing its values
# WINDOW_SIZE = config["window_size"]
# PRIORITY_ASSEMBLIES = [
#     "GCF_000001735.4",  # Arabidopsis thaliana
#     "GCF_000309985.2",  # Brassica rapa
# ]
# splits = ["train", "validation", "test"]
# EMBEDDING_WINDOW_SIZE = 100
# CHROMS = ["1", "2", "3", "4", "5"]
# NUCLEOTIDES = list("ACGT")
# models = [
#     "gonzalobenegas/gpn-brassicales",
#     #"/scratch/users/gbenegas/checkpoints/GPN_Arabidopsis_multispecies/ConvNet_ss_12k/checkpoint-12000",
# ]
# all_radius = [
#     100_000,
# ]
# all_threshold = [
#     0.1,
# ]

# # Define the tasks
# gwf.target("filter_assemblies", inputs=[config["raw_assemblies_path"]], outputs=[config["assemblies_path"]]) << f"""
# python -m gpn.ss.filter_assemblies {config['raw_assemblies_path']} {config['assemblies_path']} --keep_one_per_genus --priority_assemblies {PRIORITY_ASSEMBLIES}
# """

# # Add more tasks as needed
# # ...

# # Example of adding a task for each model
# for model in models:
#     gwf.target(f"vep_embeddings_{model}", inputs=[], outputs=[f"output/variants/all/vep_embeddings/{model}.parquet"]) << f"""
#     # Command to generate vep_embeddings for {model}
#     # ...
#     """

# # Example of adding a task for each split
# for split in splits:
#     gwf.target(f"balanced_data_{split}", inputs=[], outputs=[f"output/merged_dataset/{config['window_size']}/{config['step_size']}/{config['add_rc']}/balanced/data/{split}"]) << f"""
#     # Command to generate balanced data for {split}
#     # ...
#     """

# # Example of adding a task for each nucleotide
# for nuc in NUCLEOTIDES:
#     gwf.target(f"bed_probs_{nuc}", inputs=[], outputs=[f"output/whole_genome/bed_probs/{model}/{nuc}.bw"]) << f"""
#     # Command to generate bed_probs for {nuc}
#     # ...
#     """

# # Example of adding a task for each radius and threshold
# for radius in all_radius:
#     for threshold in all_threshold:
#         gwf.target(f"LD_{radius}_{threshold}", inputs=[], outputs=[f"output/aragwas/LD_{radius}_{threshold}.npz"]) << f"""
#         # Command to generate LD for radius {radius} and threshold {threshold}
#         # ...
#         """
