import pandas as pd
import numpy as np
import os
import glob

phylowgs_dir = "/path/to/phylowgs"
truthfiles_dir = "/path/to/groundtruth/files"
finished_samples = "/path/to/samplesfolder/finished/for/phylowgs"
samples = [dir for dir in os.listdir(finished_samples) if os.path.isdir(os.path.join(finished_samples, dir))]

output_file = os.path.join(phylowgs_dir, "statistics_plot.txt")

with open(output_file, 'w') as out_f:
    for sample in samples:
        sample_id = sample.replace("sample", "")
        print(sample_id)

        vcf_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}-noXY.truth.scoring_vcf.vcf")
        clones_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}-noXY.truth.2A.txt")
        phylowgs_res_path = os.path.join(phylowgs_dir, f"sample{sample_id}", "clonal_ssms.txt")

        try:
            with open(vcf_file_path, 'r') as vcf_file:
                header = ""
                for line in vcf_file:
                    if line.startswith("##"):
                        continue  
                    elif line.startswith("#"): 
                        header = line.strip().split("\t")
                        break 

            header.append('Status')
            df = pd.read_csv(vcf_file_path, sep='\t', comment='#', names=header, index_col=False)

            df_filtered = df[df['Status'] == True]

            clones = np.loadtxt(clones_file_path)
            df_filtered.loc[:, 'Cluster'] = clones
            df_filtered['Cluster'] = df_filtered['Cluster'].astype(int)

            #keep only cluster 1
            df_filtered = df_filtered[df_filtered['Cluster'] == 1]
            df_filtered.reset_index(drop=True, inplace=True)
            
            clonal_truth = f"{truthfiles_dir}/{sample_id}-noXY/{sample_id}_clonal_truth.txt"
            
            truth_nonclonal_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}_nonclonal_truth.txt")
            truth_nonclonal = pd.read_csv(truth_nonclonal_file_path, sep='\t', header=0, index_col=False)
            
            #df_filtered.to_csv(clonal_truth, sep='\t', index=False)

            #read Phylowgs results
            phylowgs_res = pd.read_csv(phylowgs_res_path, sep='\t', header=0, index_col=False)
            phylowgs_res[['#CHROM', 'POS']] = phylowgs_res['gene'].str.split('_', expand=True)
            phylowgs_res['#CHROM'] = phylowgs_res['#CHROM'].astype(int)
            phylowgs_res['POS'] = phylowgs_res['POS'].astype(int)
            phylowgs_res = phylowgs_res.sort_values(by=['#CHROM', 'POS'])
            
            phylowgs_mut = f"{phylowgs_dir}/sample{sample_id}/{sample_id}_clonalmutations.txt"
            phylowgs_res.to_csv(phylowgs_mut, sep='\t', index=False)

            #find common mutations
            common_mut = phylowgs_res.merge(df_filtered, on=['#CHROM', 'POS'])
            stat = len(common_mut) / len(df_filtered) * 100
            false_positives = phylowgs_res.merge(df_filtered, on=['#CHROM', 'POS'], how='left', indicator=True)
            false_positives = false_positives[false_positives['_merge'] == 'left_only']
            
            tn_merge = truth_nonclonal.merge(phylowgs_res, on=['#CHROM', 'POS'], how='left', indicator=True)
            true_negatives = tn_merge[tn_merge['_merge'] == 'left_only']
        
            fpr = (len(false_positives) / (len(false_positives) + len(true_negatives))) * 100
            
            clonal_mut_number = len(df_filtered)

            out_f.write(f"Percent of identified mutations by Phylowgs in {sample_id}: {stat:.2f}%\n")
            out_f.write(f"Percent of false positives in {sample_id}: {fpr:.2f}%\n")
            out_f.write(f"True clonal mutations in {sample_id}: {clonal_mut_number}\n")
            print(f"Processed {sample_id}: {stat:.2f}%")

        except Exception as e:
            print(f"Error processing {sample_id}: {e}")
            out_f.write(f"Error processing {sample_id}\n")

