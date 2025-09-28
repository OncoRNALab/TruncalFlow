import os
import pandas as pd

base_dir = "/path/to/fastclone"

samples = [f"P{i}" for i in range(1, 26)] + [f"S{i}" for i in range(2, 11)] + [f"T{i}" for i in range(0, 17)]

results = []

for sample in samples:
    file_path = os.path.join(base_dir, sample, "fastclone_result", "subclones.csv")
    
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        num_clusters = df.shape[0]  
        results.append((sample, num_clusters))
    else:
        print(f"File not found: {file_path}")

output_file = "fastclone_clusters_per_sample.txt"
with open(output_file, "w") as f:
    f.write("Sample\tNumClusters\n")
    for sample, num_clusters in results:
        f.write(f"{sample}\t{num_clusters}\n")
