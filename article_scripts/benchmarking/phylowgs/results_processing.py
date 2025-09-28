import json
import os

finished_samples = "/path/to/samplesfolder/finished/for/phylowgs"
samples = [dir for dir in os.listdir(finished_samples) if os.path.isdir(os.path.join(finished_samples, dir))]
print(samples)
json_numbers = []
base_json_path = "/path/to/phylowgs/witness/data"

for sample in samples:
    trees_sum = open(f"{base_json_path}/{sample}/{sample}.summ.json")

    data = json.load(trees_sum)
    
    best_tree = ''
    best_score = -1000000000
    
    for tree in data['trees']:
        score = data['trees'][tree]['llh']
        if score > best_score:
            best_score = score
            best_tree = tree
    
    json_numbers.append(best_tree)

    trees_sum.close()


base_path = "/path/to/phylowgs"

for sample, tree_number in zip(samples, json_numbers):
    json_file = f"{base_json_path}/{sample}/{tree_number}.json"
    ssm_file = f"{base_path}/{sample}/prep_files/ssm_data.txt"
    output_file = f"{base_path}/{sample}/clonal_ssms.txt"

    with open(json_file, "r") as f:
        data = json.load(f)

    with open(ssm_file, "r") as f:
        ssm_data = f.readlines()

    header = ssm_data[0]  # first line is the header
    
    ssm_data = ssm_data[1:]

    #get the SSMs assigned to cluster 2 = clonal
    clonal_ssms = data["mut_assignments"]["1"]["ssms"]

    subset_ssms = [header] + [ssm_data[int(ssm[1:])] for ssm in clonal_ssms]

    with open(output_file, "w") as f:
        f.writelines(subset_ssms)
