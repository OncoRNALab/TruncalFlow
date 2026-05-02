import os
import glob
import re
import sys
import subprocess

def submit_clustering_jobs(output_dir, time, executor='local', cpus=4, mem='16G'):
    """
    Submits QuantumClone clustering jobs for each sample in output_dir.
    
    Args:
        output_dir (str): Path to directory containing sample subdirectories.
        time (str): PBS walltime (ignored in local mode).
        executor (str): 'local' or 'pbs'.
        cpus (int): CPU cores.
        mem (str): Memory limit.
    """
    output_dir = os.path.abspath(output_dir)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    clustering_script = os.path.join(script_dir, "clustering.R")
    sif_path = os.path.join(script_dir, "thesis_2025_marina_latest.sif")

    sample_folders = [d for d in glob.glob(os.path.join(output_dir, "*")) if os.path.isdir(d)]
    submitted_job_ids = []

    for sample in sample_folders:
        sample_name = os.path.basename(sample)
        jobs_dir = os.path.join(sample, "pbs_jobs")
        os.makedirs(jobs_dir, exist_ok=True)

        snv_file = glob.glob(os.path.join(sample, "*_SNVlist.txt"))
        if not snv_file:
            print(f"Skipping {sample_name} - No SNVlist file found.")
            continue

        freec_files = glob.glob(os.path.join(sample, "*_freec.txt"))
        if not freec_files:
            print(f"[WARN] No CNV file found for {sample_name}. Proceeding without CNV integration.")
        else:
            freec_file = freec_files[0]

        if executor == 'local':
            #LOCAL MODE: run directly
            out_file = os.path.join(jobs_dir, f"RunClustering_{sample_name}.out")
            err_file = os.path.join(jobs_dir, f"RunClustering_{sample_name}.err")

            cmd = [
                'apptainer', 'exec',
                f'--cpus={cpus}',
                f'--memory={mem}',
                '--cleanenv',
                '--no-home',
                '--env', 'R_LIBS_USER=/usr/local/lib/R/library',
                sif_path,
                'Rscript', clustering_script, output_dir, sample_name
            ]
            
            print(f"[INFO] Running local clustering for {sample_name}: {' '.join(cmd)}")
            result = subprocess.run(cmd, cwd=sample, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            with open(out_file, 'w') as f:
                f.write(result.stdout)
            with open(err_file, 'w') as f:
                f.write(result.stderr)

            if result.returncode == 0:
                print(f"[INFO] Local clustering complete for {sample_name}")
            else:
                print(f"[ERROR] Local clustering failed for {sample_name}")
                continue

        elif executor == 'pbs':
            #PBS MODE: generate job script
            job_script_path = os.path.join(jobs_dir, f"job_{sample_name}.sh")
            
            with open(job_script_path, "w") as fh:
                fh.write(f"""#!/bin/bash
#PBS -N QC_{sample_name}
#PBS -l nodes=1:ppn={cpus}
#PBS -l walltime={time}
#PBS -l mem={mem}
#PBS -o {jobs_dir}/RunClustering_{sample_name}.out
#PBS -e {jobs_dir}/RunClustering_{sample_name}.err
#PBS -d {sample}

apptainer exec -e --no-home --env R_LIBS_USER=/usr/local/lib/R/library {sif_path} Rscript {clustering_script} "{output_dir}" "{sample_name}"
""")

            submit_cmd = f"qsub {job_script_path}"
            result = subprocess.run(submit_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"[ERROR] Failed to submit PBS job for {sample_name}: {result.stderr}")
                continue
            else:
                job_id = result.stdout.strip().split('.')[0]
                print(f"[INFO] Submitted PBS job {job_id} for {sample_name}")
                submitted_job_ids.append(job_id)

        else:
            raise ValueError(f"Unknown executor: {executor}")

    return submitted_job_ids
            

def check_clustering_errors(output_dir):
    sample_dirs = [
        d for d in os.listdir(output_dir)
        if os.path.isdir(os.path.join(output_dir, d))
    ]
    failed_samples = []

    for sample in sample_dirs:
        log_dir = os.path.join(output_dir, sample, "pbs_jobs")
        err_file = os.path.join(log_dir, f"RunClustering_{sample}.err")

        if not os.path.exists(err_file):
            continue

        with open(err_file, "r") as f:
            content = f.read().lower()

        if re.search(r"\b(error|oom|killed|exception|traceback)\b", content):
            print(f"[ERROR] Clustering failed for sample {sample}. Check: {err_file}")
            failed_samples.append(sample)

    return failed_samples


