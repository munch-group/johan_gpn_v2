# %% [markdown]
# ---
# title: GWF workflow
# execute:
#   eval: false
# ---

# %% [markdown]
"""
## Imports and utility functions
"""

# %%
import os
from pathlib import Path
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
import glob
import yaml
import pandas as pd
import re

# %% [markdown]
"""
Instantiate the workflow with the name of the project folder:
"""

# %%
# instantiate the workflow
gwf = Workflow(defaults={'account': 'baboons'})


config = yaml.safe_load(open("workflow.yaml"))

assemblies = pd.read_csv(config['assemblies_path'], sep='\t')
assemblies["Assembly Name"] = assemblies["Assembly Name"].str.replace(" ", "_")
assemblies.set_index("Assembly Accession", inplace=True)
assemblies["genome_path"] = [f"steps/tmp/{i}/ncbi_dataset/data/{i}/{i}_{n}_genomic.fna" for i, n in zip(assemblies.index, assemblies["Assembly Name"])]
assemblies["annotation_path"] = [f"steps/tmp/{i}/ncbi_dataset/data/{i}/genomic.gff" for i in assemblies.index]

# %% [markdown]
"""
Utility functions:
"""

# %%
# utility function
def modify_path(path, **kwargs):
    """
    Utility function for modifying file paths substituting
    the directory (dir), base name (base), or file suffix (suffix).
    """
    for key in ['dir', 'base', 'suffix']:
        kwargs.setdefault(key, None)
    assert len(kwargs) == 3

    par, name = os.path.split(path)
    name_no_suffix, suf = os.path.splitext(name)
    if type(kwargs['suffix']) is str:
        suf = kwargs['suffix']
    if kwargs['dir'] is not None:
        par = kwargs['dir']
    if kwargs['base'] is not None:
        name_no_suffix = kwargs['base']

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(kwargs['suffix']) is tuple:
        assert len(kwargs['suffix']) == 2
        new_path, nsubs = re.subn(r'{}$'.format(kwargs['suffix'][0]), kwargs['suffix'][1], new_path)
        assert nsubs == 1, nsubs
    return new_path

# %% [markdown]
"""
## Template functions:
"""
# %%

# task template function
def download_genome(assembly): 
    """
    Formats names to uppercase.
    """

    tmp_dir=f"steps/tmp/{assembly}"
    genome_path = assemblies.loc[assembly, "genome_path"]
    annotation_path = assemblies.loc[assembly, "annotation_path"]

    inputs = [config['assemblies_path']]
    genome_file = f"steps/genome/{assembly}.fa.gz"
    annotation_file = f"steps/annotation/{assembly}.gff.gz"
    outputs = [genome_file, annotation_file]
    # resource specification
    options = {'memory': '8g', 'walltime': '02:00:00'} 

    spec = f"""
        orig=$(pwd)
        mkdir -p steps/genome && 
        mkdir -p steps/annotation && 
        mkdir -p {tmp_dir} && 
        cd {tmp_dir} && 
        datasets download genome accession {assembly} --include genome,gff3 &&
        unzip -o ncbi_dataset.zip && 
        cd $orig && 
        gzip -c {genome_path} > {genome_file} && 
        gzip -c {annotation_path} > {annotation_file} && rm -r {tmp_dir}
    """
    # return target
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_all_intervals(assembly):
    inputs = [f"steps/genome/{assembly}.fa.gz"]
    outputs = [f"steps/intervals/{assembly}/all.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&
    python scripts/make_all_intervals.py {inputs[0]} {outputs[0]} {config['window_size']}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_defined_intervals(assembly):
    inputs = [f"steps/genome/{assembly}.fa.gz"]
    outputs = [f"steps/intervals/{assembly}/defined.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&
    python scripts/make_defined_intervals.py {inputs[0]} {outputs[0]} {config['window_size']}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_annotation_intervals(assembly, feature):
    inputs = [f"steps/intervals/{assembly}/defined.parquet",
              f"steps/genome/{assembly}.fa.gz"]
    outputs = [f"steps/intervals/{assembly}/annotation_{feature}.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    include_flank = config.get("annotation_features_include_flank", config['window_size'] // 2)
    add_jiter = config.get("annotation_features_add_jitter", 100)
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&
    python scripts/make_annotation_intervals.py {inputs[0]} {inputs[1]} {outputs[0]} \
        {config['window_size']} {feature} {include_flank} {add_jiter}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_balanced_v1_intervals(assembly):
    inputs = [f"steps/intervals/{assembly}/defined.parquet",
              f"steps/annotation/{assembly}.gff.gz"]
    outputs = [f"steps/intervals/{assembly}/balanced_v1.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    promoter_upstream = config.get("promoter_upstream", 1000)
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&
    python scripts/make_defined_intervals.py {inputs[0]} {outputs[0]} \
        {config['window_size']} {promoter_upstream}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_dataset_assembly(assembly):
    splits = ["train", "validation", "test"]
    inputs = [f"results/intervals/{assembly}/{config['target_intervals']}.parquet",
              f"results/genome/{assembly}.fa.gz"]
    outputs = [f"results/dataset_assembly/{assembly}/{split}.parquet" for split in splits]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&    
    python scripts/make_dataset_assembly.py {' '.join(inputs)} {' '.join(outputs)} \
        {config['split_proportion']} {config['window_size']} {config['step_size']} {config['add_rc']} \
        {config['whitelist_validation_chroms']} {config['whitelist_test_chroms']}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


download_targets = gwf.map(download_genome, assemblies.index)

if config['target_intervals'] == 'all':
    interval_targets = gwf.map(make_all_intervals, assemblies.index)
elif config['target_intervals'] == 'defined':
    interval_targets = gwf.map(make_defined_intervals, assemblies.index)
elif config['target_intervals'].startswith('annotation'):
    feature = config['target_intervals'].replace('annotation_', '')
    interval_targets = gwf.map(make_annotation_intervals, assemblies.index, feature)
elif config['target_intervals'] == 'balanced_v1':
    interval_targets = gwf.map(make_balanced_v1_intervals, assemblies.index)
else:
    assert 0

datasets = gwf.map(make_dataset_assembly, assemblies.index)


# from gpn.data import (
#     Genome, load_table, get_balanced_intervals, filter_length,
#     filter_annotation_features,
# )


# # before uploading to HF Hub, remove data/split/.snakemake_timestamp files
# rule merge_datasets:
#     input:
#         expand("results/dataset_assembly/{assembly}/{{split}}.parquet", assembly=assemblies.index),
#     output:
#         directory("results/dataset/data/{split}"),
#     threads: workflow.cores
#     run:
#         intervals = pd.concat(
#             tqdm((pd.read_parquet(path) for path in input), total=len(input)),
#             ignore_index=True,
#         ).sample(frac=1, random_state=42)
#         print(intervals)

#         if config.get("subsample_to_target", False) and wildcards.split == "train":
#             n_target = (intervals.assembly==config["target_assembly"]).sum()
#             intervals = intervals.groupby("assembly").sample(
#                 n=n_target, random_state=42
#             ).sample(frac=1, random_state=42)
#             print(wildcards.split, intervals.assembly.value_counts())
#             print(intervals)

#         n_shards = math.ceil(len(intervals) / config["samples_per_file"])
#         assert n_shards < 10000
#         os.makedirs(output[0])
#         for i in tqdm(range(n_shards)):
#             path = Path(output[0]) / f"shard_{i:05}.jsonl.zst"
#             intervals.iloc[i::n_shards].to_json(
#                 path, orient="records", lines=True,
#                 compression={'method': 'zstd', 'threads': -1}
#             )



# # %%

# %%
