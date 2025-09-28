import pandas as pd
import glob
import os
import csv

battenberg_file = glob.glob('*battenberg.txt')[0]

def get_sample_id(directory):
    parts = os.path.basename(directory).split('-')
    if len(parts) > 1:
        return parts[0]
    return None

current_dir = os.getcwd()
sample_id = get_sample_id(current_dir)

header = ["CHROM", "START", "END"]


loh_regions = []
with open(battenberg_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if int(row['nMin1_A']) == 0:  #check for LOH
            loh_regions.append([row['chr'], row['startpos'], row['endpos']])

loh_file = pd.DataFrame(loh_regions, columns=header)

base_output_dir = "/path/to/sciclonefolder"
output_dir = os.path.join(base_output_dir, f"sample{sample_id}")

output_filename = os.path.join(output_dir, f"{sample_id}_loh.txt")

loh_file.to_csv(output_filename, sep='\t', index=False)

