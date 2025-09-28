import pandas as pd
import numpy as np
import os
import glob


quantum_dir = "/path/to/quantumclonefolder"
truthfiles_dir = "/path/to/groundtruth/files"
finished_samples = "/path/to/quantumclonefolder"
samples = [dir for dir in os.listdir(finished_samples) if os.path.isdir(os.path.join(finished_samples, dir))]
output_file = "/path/to/quantumclonefolder/statistics_quantumclone_allsamples_plot.txt"

with open(output_file, 'w') as out_f:
    for sample in samples:
        sample_id = sample.replace("sample", "")
     
        truth_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}_clonal_truth.txt")
        truth_nonclonal_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}_nonclonal_truth.txt")
        
        filtered_path = os.path.join(quantum_dir, f"sample{sample_id}", "filtered.csv")
        cluster_path = os.path.join(quantum_dir, f"sample{sample_id}", "clustering.csv")
        centers_path = os.path.join(quantum_dir, f"sample{sample_id}", "centers.csv")
        
        filtered = pd.read_csv(filtered_path, sep=',', header=0, index_col=0)
        clusters = pd.read_csv(cluster_path, sep=',', header=0, index_col=0)
        centers = pd.read_csv(centers_path, sep=',', header=0, index_col=0)
        
        clonal_cluster = int(centers["X..i.."].idxmax())
        
        filtered = filtered.sort_values(by="id")
        
        filtered['Cluster'] = clusters['Number']
        
        quantumclone_clonal = filtered[filtered['Cluster'] == clonal_cluster]
        
        quantumclone_clonal = quantumclone_clonal.rename(columns={'Chr': '#CHROM', 'Start': 'POS'})
        
        quantumclone_mut = f"{quantum_dir}/sample{sample_id}/{sample_id}_clonalmutations.txt"
        quantumclone_clonal.to_csv(quantumclone_mut, sep='\t', index=False)

        truth = pd.read_csv(truth_file_path, sep = '\t', index_col=False)
        truth_nonclonal = pd.read_csv(truth_nonclonal_file_path, sep='\t', header=0, index_col=False)
        
        common_mut = quantumclone_clonal.merge(truth, on=['#CHROM', 'POS'])
        false_positives = quantumclone_clonal.merge(truth, on=['#CHROM', 'POS'], how='left', indicator=True)
        false_positives = false_positives[false_positives['_merge'] == 'left_only']
        tn_merge = truth_nonclonal.merge(quantumclone_clonal, on=['#CHROM', 'POS'], how='left', indicator=True)
        true_negatives = tn_merge[tn_merge['_merge'] == 'left_only']
    
        fpr = (len(false_positives) / (len(false_positives) + len(true_negatives))) * 100
        
        stat = len(common_mut) / len(truth) * 100
        
        clonal_mut_number = len(truth)

        out_f.write(f"Percent of identified mutations by QuantumClone in {sample_id}: {stat:.2f}%\n")
        out_f.write(f"Percent of false positives in {sample_id}: {fpr:.2f}%\n")
        out_f.write(f"True clonal mutations in {sample_id}: {clonal_mut_number}\n")
