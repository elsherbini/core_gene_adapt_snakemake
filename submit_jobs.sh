#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH --partition=medium
#SBATCH --job-name=submit_jobs

mkdir -p ./slurm

snakemake --use-conda --conda-prefix=~/snakemake_conda_prefix --cluster-config cluster.yaml \
  --cluster "sbatch --output=slurm/{rulename}.{jobid} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem} --time={cluster.time}" \
  --jobname {rulename}.{jobid} --jobs 500 --keep-going --groups rename_fastas_ffn=group_b rename_fastas_faa=group_c --group-components group_b=10 group_c=10 --rerun-incomplete --latency-wait 20