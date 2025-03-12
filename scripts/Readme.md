# Workflow Explanation (workflow.py)

This workflow script is designed to process genomic data using the GWF workflow management system. The script defines several functions to download genomes, create intervals, and assemble datasets. Below is an explanation of the key components and functions in the script:

## Key Components

1. **Imports and Configuration:**
   - The script imports necessary libraries and modules, including `os`, `pathlib`, `gwf`, `yaml`, `pandas`, and others.
   - The workflow is instantiated with the name of the project folder.
   - Configuration settings are loaded from a YAML file (`workflow.yaml`).
   - The base directory for data storage is defined.
   - Genomic assemblies are loaded and processed using pandas.

2. **Utility Functions:**
   - `modify_path`: A utility function for modifying file paths by substituting the directory, base name, or file suffix.

3. **Task Template Functions:**
   - `download_genome`: Downloads genome and annotation files for a given assembly using the NCBI datasets command-line tool.
   - `make_all_intervals`: Creates intervals for all positions in the genome.
   - `make_defined_intervals`: Creates intervals for positions with defined nucleotides (not N).
   - `make_annotation_intervals`: Creates intervals for specific annotation features (e.g., CDS, exon).
   - `make_balanced_v1_intervals`: Creates balanced intervals based on a specific recipe.
   - `make_dataset_assembly`: Assembles datasets for training, validation, and testing.

4. **Workflow Execution:**
   - The script maps the `download_genome` function to all assemblies.
   - Based on the `target_intervals` configuration, the appropriate interval creation function is mapped to all assemblies.
   - The script prints the value of `target_intervals` for debugging purposes.

## Example Usage

To run the workflow, you need to ensure that the configuration file (`workflow.yaml`) is correctly set up. You can then use the following commands to clean and run the workflow:

```bash
gwf clean --all
gwf run