import pandas as pd
import glob
import os
import csv

vcf_file = glob.glob('*mutect.vcf')[0]
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
#output['Genotype'] = 'AB'

output['Depth'] = output['Depth'].astype(int)
output['Alt'] = output['Alt'].astype(int)


base_output_dir = "/path/to/quantumclonecnvfolder"
output_dir = os.path.join(base_output_dir, f"sample{sample_id}")

output_filename = os.path.join(output_dir, f"{sample_id}_quantumclone.txt")

output.to_csv(output_filename, sep='\t', index=False)


def convert_battenberg_to_freec(battenberg_file):
    
    header = ["Chromosome", "Start", "Ratio", "MedianRatio", "CopyNumber", "BAF", "EstimatedBAF", "Genotype", "UncertaintyOfGT"]

    
    with open(battenberg_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        freec = []
        for row in reader:
            
            ratio = 2 ** float(row['LogR'])
            
            copy_number = 0
            
            corr_baf = abs(round(float(row['BAF']), 2) - 0.5)
            
            if int(row['nMaj1_A']) == 0 and int(row['nMin1_A']) == 0:
                genotype = "0"  # Homozygous deletion in major clone
            elif int(row['nMaj1_A']) == 1 and int(row['nMin1_A']) == 0:
                genotype = "A"  # Hemizygous deletion
            else:
                genotype = "A" * int(row['nMaj1_A']) + "B" * int(row['nMin1_A'])
            
            if row['frac2_A'] != 'NA' and float(row['frac2_A']) > 0:
                cn = (int(row['nMaj1_A']) + int(row['nMin1_A'])) * float(row['frac1_A']) + (int(row['nMaj2_A']) + int(row['nMin2_A'])) * float(row['frac2_A'])
                if round(cn) == 0:
                    genotype = "0"
                elif round(cn) == 1:
                    genotype = "A"
                elif round(cn) == 2:
                    if float(row['frac1_A']) > float(row['frac2_A']):
                        if int(row['nMaj1_A']) > int(row['nMin1_A']):
                            genotype = "AA"
                        else:
                            genotype = "AB"
                    else:
                        if int(row['nMaj2_A']) > int(row['nMin2_A']):
                            genotype = "AA"
                        else:
                            genotype = "AB"

            state1_cn = int(row['nMaj1_A']) + int(row['nMin1_A'])
            copy_number += state1_cn * float(row['frac1_A'])
            
            if row['nMaj2_A'] != 'NA' and row['nMaj2_A'] != 'NA' and row['frac2_A'] != 'NA':
                state2_cn = int(row['nMaj2_A']) + int(row['nMin2_A'])
                copy_number += state2_cn * float(row['frac2_A'])
                
            copy_number = round(copy_number)
            
            if row['SDfrac_A'] != "NA":
                uncertainty = round(100 * float(row['SDfrac_A']))
            else:
                uncertainty = 0
        
            freec.append([row['chr'], row['startpos'], ratio, ratio, copy_number, corr_baf, round(float(row['BAF']), 2), genotype, uncertainty])
            
            freec_file = pd.DataFrame(freec, columns=header)
    
    return freec_file
    
freec = convert_battenberg_to_freec(battenberg_file)

freec_filename = os.path.join(output_dir, f"{sample_id}_freec.txt")

freec.to_csv(freec_filename, sep='\t', index=False)

