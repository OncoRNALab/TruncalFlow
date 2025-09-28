import pandas as pd
import os

def parse_vcf(vcf_file):
    mutations = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("##"):
                continue  
            elif line.startswith("#"): 
                continue
            else:
                parts = line.strip().split("\t")
                chrom, pos, _, ref, _, _, _, _, fmt, _, tumor = parts
                
                # Extract ref and alt counts from AD field
                fmt_fields = fmt.split(":")
                tumor_fields = tumor.split(":")
                
                if "AD" in fmt_fields:
                    ad_index = fmt_fields.index("AD")
                    ref_counts, alt_counts = map(int, tumor_fields[ad_index].split(","))
                else:
                    continue  #skip if AD field is missing
                
                mutation_id = f"{chrom}:{pos}:{ref}"
                mutations.append([mutation_id, chrom, int(pos), ref_counts, alt_counts])
    
    return pd.DataFrame(mutations, columns=["mutation_id", "chr", "pos", "ref_counts", "alt_counts"])


#reading Battenberg CNV File
def parse_battenberg(cnv_file):
    cnv_data = pd.read_csv(cnv_file, sep="\t")
    
    cnv_data = cnv_data[["chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A"]]
    cnv_data["normal_cn"] = 2  # Autosomal default

    return cnv_data


#merging VCF and CNV Data
def merge_data(vcf_df, cnv_df, sample_id):
    result = []
    
    cnv_df["chr"] = cnv_df["chr"].astype(str)
    
    for _, row in vcf_df.iterrows():
        row["chr"] = str(row["chr"])
        chr_match = cnv_df[cnv_df["chr"] == row["chr"]]
        match = chr_match[(chr_match["startpos"] <= row["pos"]) & (chr_match["endpos"] >= row["pos"])]
        
        if not match.empty:
            match_row = match.iloc[0]
            result.append([
                row["mutation_id"], sample_id, row["ref_counts"], row["alt_counts"],
                int(match_row["normal_cn"]), int(match_row["nMaj1_A"]), int(match_row["nMin1_A"]),
                match_row["frac1_A"]
            ])
    
    return pd.DataFrame(result, columns=["mutation_id", "sample_id", "ref_counts", "alt_counts",
                                         "normal_cn", "major_cn", "minor_cn", "tumour_content"])


samples = [f"P{i}" for i in range(1, 26)] + \
          [f"S{i}" for i in range(2, 11)] + \
          [f"T{i}" for i in range(0, 17)]

vcf_dir = "/path/to/groundtruth/files"
output_dir = "/path/to/pyclonefolder"

for sample in samples:
    vcf_file = os.path.join(vcf_dir, f"{sample}-noXY", f"{sample}-noXY.mutect.vcf")
    cnv_file = os.path.join(vcf_dir, f"{sample}-noXY", f"{sample}-noXY.battenberg.txt")
    
    if os.path.exists(vcf_file) and os.path.exists(cnv_file):
        vcf_data = parse_vcf(vcf_file)
        cnv_data = parse_battenberg(cnv_file)
        final_df = merge_data(vcf_data, cnv_data, sample)
        
        output_file = os.path.join(output_dir, f"sample{sample}", f"{sample}_pyclone_input.tsv")
        final_df.to_csv(output_file, sep="\t", index=False)
    else:
        print(f"Skipping {sample}, missing files!")
