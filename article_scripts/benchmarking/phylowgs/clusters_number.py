import json
import os

finished_samples = "/path/to/samplesfolder/finished/for/phylowgs"
samples = [dir for dir in os.listdir(finished_samples) if os.path.isdir(os.path.join(finished_samples, dir))]
json_numbers = []
base_json_path = "/path/to/phylowgs/witness/data"

base_path = "/path/to/phylowgs"
cluster_count_output = "phylowgs_clusters_per_sample.txt"
cluster_counts = []

for sample in samples:
    trees_sum = open(f"{base_json_path}/{sample}/{sample}.summ.json")

    data = json.load(trees_sum)
    
    best_tree = ''
    best_score = float('-inf')
    
    for tree in data['trees']:
        score = data['trees'][tree]['llh']
        if score > best_score:
            best_score = score
            best_tree = tree
    
    json_numbers.append(best_tree)

    trees_sum.close()

for sample, tree_number in zip(samples, json_numbers):
    json_file = f"{base_json_path}/{sample}/{tree_number}.json"
    ssm_file = f"{base_path}/{sample}/prep_files/ssm_data.txt"
    output_file = f"{base_path}/{sample}/clonal_ssms.txt"

    with open(json_file, "r") as f:
        data = json.load(f)
        
    num_clusters = len(data["mut_assignments"])
    cluster_counts.append((sample, num_clusters))

    with open(ssm_file, "r") as f:
        ssm_data = f.readlines()

    header = ssm_data[0]  # first line is the header
    
    ssm_data = ssm_data[1:]

    #get the SSMs assigned to cluster 2 = clonal
    clonal_ssms = data["mut_assignments"]["1"]["ssms"]

    subset_ssms = [header] + [ssm_data[int(ssm[1:])] for ssm in clonal_ssms]

    with open(output_file, "w") as f:
        f.writelines(subset_ssms)
    
        
with open(cluster_count_output, "w") as f:
    f.write("Sample\tNumClusters\n")
    for sample, count in cluster_counts:
        f.write(f"{sample}\t{count}\n")
