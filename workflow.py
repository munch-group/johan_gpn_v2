# %% [markdown]
# ---
# title: Making a dataset for GPN 
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
import time

# %% [markdown]
"""
Instantiate the workflow with the name of the project folder:
"""

# %%
# instantiate the workflow
#gwf = Workflow(defaults={'account': 'baboons'})
gwf = Workflow(defaults={'account': 'TopicsInBioinformatics'})


config = yaml.safe_load(open("scripts/generate_dataset/workflow.yaml"))

# Define the base directory
base_dir = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn"



## loading data and preparatoin (pandas)
assemblies = pd.read_csv(config['assemblies_path'], sep='\t')
assemblies["Assembly Name"] = assemblies["Assembly Name"].str.replace(" ", "_")
assemblies.set_index("Assembly Accession", inplace=True)
assemblies["genome_path"] = [f"{base_dir}/steps/tmp/{i}/ncbi_dataset/data/{i}/{i}_{n}_genomic.fna" for i, n in zip(assemblies.index, assemblies["Assembly Name"])]
assemblies["annotation_path"] = [f"{base_dir}/steps/tmp/{i}/ncbi_dataset/data/{i}/genomic.gff" for i in assemblies.index]

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
    #print(f"Downloading genome for assembly: {assembly}")

    tmp_dir = f"{base_dir}/steps/tmp/{assembly}"
    genome_path = assemblies.loc[assembly, "genome_path"]
    annotation_path = assemblies.loc[assembly, "annotation_path"]

    inputs = [config['assemblies_path']]
    genome_file = f"{base_dir}/steps/genome/{assembly}.fa.gz"
    annotation_file = f"{base_dir}/steps/annotation/{assembly}.gff.gz"
    outputs = [genome_file, annotation_file]
    options = {'memory': '8g', 'walltime': '02:00:00'}

    spec = f"""
        orig=$(pwd)
        mkdir -p {base_dir}/steps/genome && 
        mkdir -p {base_dir}/steps/annotation && 
        mkdir -p {tmp_dir} && 
        cd {tmp_dir} && 
        datasets download genome accession {assembly} --include genome,gff3 &&
        unzip -o ncbi_dataset.zip && 
        cd $orig && 
        gzip -c {genome_path} > {genome_file} && 
        gzip -c {annotation_path} > {annotation_file} && rm -r {tmp_dir}
    """

    #print(f"Spec for {assembly}: {spec}")
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_all_intervals(assembly):
    inputs = [f"{base_dir}/steps/genome/{assembly}.fa.gz"]
    outputs = [f"{base_dir}/steps/intervals/{assembly}/all.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    spec = f"""
    mkdir -p {base_dir}/steps/intervals/{assembly} &&
    python scripts/generate_dataset/make_all_intervals.py {inputs[0]} {outputs[0]} {config['window_size']}
    """
    #print(f"Spec for make_all_intervals {assembly}: {spec}")
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_defined_intervals(assembly): ## does not show up in gwf
    inputs = [f"{base_dir}/steps/genome/{assembly}.fa.gz"]
    outputs = [f"{base_dir}/steps/intervals/{assembly}/defined.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&
    python scripts/generate_dataset/make_defined_intervals.py {inputs[0]} {outputs[0]} {config['window_size']}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_annotation_intervals(assembly, feature):
    inputs = [f"{base_dir}/steps/intervals/{assembly}/defined.parquet",
              f"{base_dir}/steps/genome/{assembly}.fa.gz"]
    outputs = [f"{base_dir}/steps/intervals/{assembly}/annotation_{feature}.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    include_flank = config.get("annotation_features_include_flank", config['window_size'] // 2)
    add_jiter = config.get("annotation_features_add_jitter", 100)
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&
    python scripts/generate_dataset/make_annotation_intervals.py {inputs[0]} {inputs[1]} {outputs[0]} \
        {config['window_size']} {feature} {include_flank} {add_jiter}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_balanced_v1_intervals(assembly):  ### maybe not that inmpotent becuase it was the date created for the paper 
    inputs = [f"{base_dir}/steps/intervals/{assembly}/defined.parquet",
              f"{base_dir}/steps/annotation/{assembly}.gff.gz"]
    outputs = [f"{base_dir}/steps/intervals/{assembly}/balanced_v1.parquet"]
    options = {'memory': '8g', 'walltime': '02:00:00'} 
    promoter_upstream = config.get("promoter_upstream", 1000)
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&
    python scripts/generate_dataset/make_defined_intervals.py {inputs[0]} {outputs[0]} \
        {config['window_size']} {promoter_upstream}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_dataset_assembly(assembly):
    splits = ["train", "validation", "test"]
    inputs = [f"{base_dir}/steps/intervals/{assembly}/{config['target_intervals']}.parquet",      ### it ask for data in a folder that does not exist
              f"{base_dir}/steps/genome/{assembly}.fa.gz"]                                        ### it ask for data in a folder that does not exist
    outputs = [f"{base_dir}/steps/dataset_assembly/{assembly}/{split}.parquet" for split in splits]
    options = {'memory': '24g', 'walltime': '02:00:00'} 
    spec = f"""
    mkdir -p steps/intervals/{assembly} &&    
    python scripts/generategenerate_dataset_data_set/make_dataset_assembly.py {' '.join(inputs)} {' '.join(outputs)} \
        {config['split_proportion']} {config['window_size']} {config['step_size']} {config['add_rc']} \
        {config['whitelist_validation_chroms']} {config['whitelist_test_chroms']}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


download_targets = gwf.map(download_genome, assemblies.index)







from gpn.data import (
    Genome, load_table, get_balanced_intervals, filter_length,
    filter_annotation_features,
)

def merge_datasets(assembly):
    splits = ["train", "validation", "test"]
    inputs = [f"{base_dir}/steps/dataset_assembly/{assembly}/{split}.parquet" for split in splits]
    output_dir = f"{base_dir}/steps/dataset/data/{assembly}"
    options = {'memory': '32g', 'walltime': '02:00:00'} 
    spec = f"""
    mkdir -p {output_dir} &&    
    python /faststorage/project/johan_gpn/people/johanulsrup/johan_gpn/scripts/generate_dataset/make_merge_datasets.py {' '.join(inputs)} {output_dir}
    """
    return AnonymousTarget(inputs=inputs, outputs=[output_dir], options=options, spec=spec)


# %% [markdown]
"""
## Logic for target intervals defind in the yaml file:
# Intervals from fasta file used for training:
# - "all": all positions
# - "defined": positions with defined nucleotides (not N)
# - "annotation_{feature}": only <feature> positions from annotation, e.g. CDS, exon
# - "balanced_v1": recipe used in original paper
"""
# %%


if config['target_intervals'] == 'all':
    interval_targets = gwf.map(make_all_intervals, assemblies.index)
elif config['target_intervals'] == 'defined':
    print("Calling make_defined_intervals")
    interval_targets = gwf.map(make_defined_intervals, assemblies.index)  
elif config['target_intervals'].startswith('annotation'):
    feature = config['target_intervals'].replace('annotation_', '')
    interval_targets = gwf.map(make_annotation_intervals, assemblies.index, feature)
elif config['target_intervals'] == 'balanced_v1':
    interval_targets = gwf.map(make_balanced_v1_intervals, assemblies.index)
else:
    assert 0

datasets = gwf.map(make_dataset_assembly, assemblies.index)
merge_datasets_targets = gwf.map(merge_datasets, assemblies.index)

# # %%

# %%



