import pandas as pd
import glob
import os

vcf_file = glob.glob('*mutect.vcf')[0]

def get_sample_id(directory):
    parts = os.path.basename(directory).split('-')
    if len(parts) > 1:
        return parts[0]
    return None

current_dir = os.getcwd()
sample_id = get_sample_id(current_dir)



with open(vcf_file, 'r') as file:
    lines = []
    header = ["CHROM", "POS", "REF", "VAR",
        "NRM_REF_READS", "NRM_VAR_READS", "NRM_VAF",
        "TUM1_REF_READS", "TUM1_VAR_READS", "TUM1_VAF"]
    for line in file:
        if line.startswith('#'):
            continue
        else:
            parts = line.strip().split('\t')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]
            normal_info = parts[9].split(':')
            tumor_info = parts[10].split(':')

            nrm_ref_reads, nrm_var_reads = map(int, normal_info[1].split(','))
            nrm_vaf = float(normal_info[4])
            
            tum_ref_reads, tum_var_reads = map(int, tumor_info[1].split(','))
            tum_vaf = float(tumor_info[4])
            
            lines.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "VAR": alt,
                "NRM_REF_READS": nrm_ref_reads,
                "NRM_VAR_READS": nrm_var_reads,
                "NRM_VAF": nrm_vaf,
                "TUM1_REF_READS": tum_ref_reads,
                "TUM1_VAR_READS": tum_var_reads,
                "TUM1_VAF": tum_vaf})
        
vcf = pd.DataFrame(lines, columns=header)

base_output_dir = "/path/to/sciclonefolder"
output_dir = os.path.join(base_output_dir, f"sample{sample_id}")

output_filename = os.path.join(output_dir, f"{sample_id}_sciclone.txt")

vcf.to_csv(output_filename, sep='\t', index=False)

