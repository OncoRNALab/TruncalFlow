import pandas as pd
import numpy as np
import os
import glob


fastclone_dir = "/path/to/fastclone"
truthfiles_dir = "/path/to/groundtruth/files"
finished_samples = "/path/to/samplesfolder/finished/for/phylowgs"
samples = [dir for dir in os.listdir(finished_samples) if os.path.isdir(os.path.join(finished_samples, dir))]
output_file = os.path.join(fastclone_dir, "statistics_fastclone_partsamples_plot.txt")

with open(output_file, 'w') as out_f:
    for sample in samples:
        sample_id = sample.replace("sample", "")
        
        truth_clonal_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}_clonal_truth.txt")
        truth_nonclonal_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}_nonclonal_truth.txt")
        fastclone_res_path = os.path.join(fastclone_dir, f"{sample_id}", "fastclone_result", "scores.csv")
        fastclone_clusters_path = os.path.join(fastclone_dir, f"{sample_id}", "fastclone_result", "subclones.csv")

        fastclone_res = pd.read_csv(fastclone_res_path, sep=',', header=0, index_col=False)
        fastclone_res = fastclone_res.dropna()
        clusters = pd.read_csv(fastclone_clusters_path, sep=',', header=0, index_col=False)
        fastclone_res.rename(columns={fastclone_res.columns[0]: "mutation_id"}, inplace=True)
        clusters.rename(columns={clusters.columns[0]: "cluster"}, inplace=True)
        
        clonal_cluster = int(clusters["cluster"].idxmax())
        cluster_cols = fastclone_res.columns[1:]  #all after 'mutation_id'

        fastclone_res["Cluster"] = fastclone_res[cluster_cols].idxmax(axis=1)
        fastclone_res['Cluster'] = fastclone_res['Cluster'].astype(int)
        
        fastclone_clones = fastclone_res[fastclone_res['Cluster'] == clonal_cluster]
    
        fastclone_clones[['#CHROM', 'POS']] = fastclone_clones['mutation_id'].str.split(':', n=2, expand=True).iloc[:, :2]
        fastclone_clones['#CHROM'] = fastclone_clones['#CHROM'].astype(int)
        fastclone_clones['POS'] = fastclone_clones['POS'].astype(int)
        fastclone_clones = fastclone_clones.sort_values(by=['#CHROM', 'POS'])
        
        fastclone_mut = f"{fastclone_dir}/{sample_id}/{sample_id}_clonalmutations.txt"
        fastclone_clones.to_csv(fastclone_mut, sep='\t', index=False)
        
        truth_clonal = pd.read_csv(truth_clonal_file_path, sep='\t', header=0, index_col=False)
        truth_nonclonal = pd.read_csv(truth_nonclonal_file_path, sep='\t', header=0, index_col=False)

        #find common mutations
        common_mut = fastclone_clones.merge(truth_clonal, on=['#CHROM', 'POS'])
        stat = len(common_mut) / len(truth_clonal) * 100
        
        false_positives = fastclone_clones.merge(truth_clonal, on=['#CHROM', 'POS'], how='left', indicator=True)
        false_positives = false_positives[false_positives['_merge'] == 'left_only']
        tn_merge = truth_nonclonal.merge(fastclone_clones, on=['#CHROM', 'POS'], how='left', indicator=True)
        true_negatives = tn_merge[tn_merge['_merge'] == 'left_only']
    
        fpr = (len(false_positives) / (len(false_positives) + len(true_negatives))) * 100
        
        clonal_mut_number = len(truth_clonal)

        out_f.write(f"Percent of identified mutations by FastClone in {sample_id}: {stat:.2f}%\n")
        out_f.write(f"Percent of false positives in {sample_id}: {fpr:.2f}%\n")
        out_f.write(f"True clonal mutations in {sample_id}: {clonal_mut_number}\n")