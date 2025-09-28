import pandas as pd
import glob
import os

vcf_file = glob.glob('*mutect.vcf')[0]

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

def parse_format_field(row, field_name):
    format_fields = row['FORMAT'].split(':')
    tumor_values = row['tumor'].split(':')
    field_index = format_fields.index(field_name)
    return tumor_values[field_index]

with open(vcf_file, 'r') as file:
    lines = []
    
    for line in file:
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            header = line.strip().split('\t')
        else:
           lines.append(line.strip().split('\t'))

vcf = pd.DataFrame(lines, columns=header)

output = pd.DataFrame({
    'Sample': vcf_file.split('.')[0],
    'Chr': vcf['#CHROM'],
    'Start': vcf['POS']
})

output['Depth'] = vcf.apply(lambda row: parse_format_field(row, 'DP'), axis=1)
output['Alt'] = vcf.apply(lambda row: parse_format_field(row, 'AD').split(',')[1], axis=1)
output['Genotype'] = 'AB'

output['Depth'] = output['Depth'].astype(int)
output['Alt'] = output['Alt'].astype(int)


base_output_dir = "/path/to/quantumclonefolder"
output_dir = os.path.join(base_output_dir, f"sample{sample_id}")

output_filename = os.path.join(output_dir, f"{sample_id}_quantumclone.txt")

output.to_csv(output_filename, sep='\t', index=False)