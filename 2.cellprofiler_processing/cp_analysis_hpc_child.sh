#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --account=amc-general
#SBATCH --time=2-00:00:00
#SBATCH --output=run_CP_child-%j.out

# 1 task at 10GB RAM for the core (adjust as needed)

# activate cellprofiler environment
module load miniforge
conda init bash
conda activate fibrosis_cp_env

# input image directory passed as first argument
image_dir=$1

# run your python analysis script with the input image directory
python nbconverted/cp_analysis_hpc.py --image_dir "$image_dir"

# deactivate conda environment
conda deactivate

echo "CellProfiler analysis done for directory: $image_dir"
