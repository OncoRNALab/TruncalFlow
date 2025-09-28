import os
import pandas as pd

samples = [f"P{i}" for i in range(1, 26)] + [f"S{i}" for i in range(2, 11)] + [f"T{i}" for i in range(0, 17)]

base_path = "/path/to/quantumclonecnvfolder"
output_file = "quantumclone_cnv_clusters_per_sample.txt"
results = []

for sample in samples:
    file_path = os.path.join(base_path, f"sample{sample}", "centers.csv")
    
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        num_clusters = df.shape[0]
        results.append((sample, num_clusters))
    else:
        print(f"File not found: {file_path}")

with open(output_file, "w") as f:
    f.write("Sample\tNumClusters\n")
    for sample, count in results:
        f.write(f"{sample}\t{count}\n")
