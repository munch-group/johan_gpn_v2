import csv

# Define the output file path
fileparth = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/scripts/gpn_arabidopsis/input/assembly_list"
output_file = f"{fileparth}/genome_metadata.tsv"
# Define the header
header = [
    "Assembly Accession",
    "Assembly Name",
    "Organism Name",
    "Organism Infraspecific Names Breed",
    "Organism Infraspecific Names Strain",
    "Organism Infraspecific Names Cultivar",
    "Organism Infraspecific Names Ecotype",
    "Organism Infraspecific Names Isolate",
    "Organism Infraspecific Names Sex",
    "Annotation Name",
    "Assembly Stats Total Sequence Length",
    "Assembly Level",
    "Assembly Submission Date",
    "WGS project accession",
    "genus",
    "Priority"
]

# Define the row for Papio anubis
papio_row = [
    "GCF_000264685.3",       # Assembly Accession
    "Panu_3.0",              # Assembly Name
    "Papio anubis",          # Organism Name
    "", "", "", "", "", "",  # Empty fields for breed, strain, etc.
    "NCBI Annotation Release 101",  # Annotation Name
    "2935174152",            # Total Sequence Length
    "Scaffold",              # Assembly Level
    "2012-08-15",            # Submission Date
    "AANL01",                # WGS Project Accession
    "Papio",                 # Genus
    "1_Low"                  # Priority
]

# Write to TSV
with open(output_file, "w", newline="") as tsvfile:
    writer = csv.writer(tsvfile, delimiter="\t")
    writer.writerow(header)
    writer.writerow(papio_row)

print(f"TSV file created: {output_file}")
