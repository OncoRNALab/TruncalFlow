import pandas as pd
import glob
import os
import csv

battenberg_file = glob.glob('*battenberg.txt')[0]

def get_sample_id(directory):
    parts = os.path.basename(directory).split('-')
    print(parts)
    if len(parts) > 1:
        return parts[0]
    return None

current_dir = os.getcwd()
print(current_dir)
sample_id = get_sample_id(current_dir)
print(sample_id)

header = ["CHROM", "START", "END", "SEGMENT_MEAN"]

copy_number = []
with open(battenberg_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            copy_number.append([row['chr'], row['startpos'], row['endpos'], row['LogR']])
            

copy_number_file = pd.DataFrame(copy_number, columns=header)

base_output_dir = "/path/to/sciclonefolder"
output_dir = os.path.join(base_output_dir, f"sample{sample_id}")

output_filename = os.path.join(output_dir, f"{sample_id}_copynumber.txt")

copy_number_file.to_csv(output_filename, sep='\t', index=False)

