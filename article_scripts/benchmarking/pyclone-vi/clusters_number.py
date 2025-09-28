import os
import pandas as pd

samples = [f"P{i}" for i in range(1, 26)] + [f"S{i}" for i in range(2, 11)] + [f"T{i}" for i in range(0, 17)]

base_path = "/path/to/pyclonefolder"
output_file = "pyclone_clusters_per_sample.txt"

results = []

for sample in samples:
    file_path = os.path.join(base_path, f"sample{sample}", f"sample{sample}_result.tsv")
    
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep="\t")
        unique_clusters = df["cluster_id"].nunique()
        results.append((sample, unique_clusters))
    else:
        print(f"File not found: {file_path}")

with open(output_file, "w") as f:
    f.write("Sample\tNumClusters\n")
    for sample, num_clusters in results:
        f.write(f"{sample}\t{num_clusters}\n")
