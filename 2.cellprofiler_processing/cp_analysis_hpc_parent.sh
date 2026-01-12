#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=10:00
#SBATCH --output=cp_parent-%j.out

# activate cellprofiler environment
module load miniforge
conda init bash
conda activate fibrosis_cp_env

# convert all notebooks to python scripts (if any exist)
jupyter nbconvert --to=script --FilesWriter.build_directory=nbconverted/ *.ipynb

# define the batch variable
batch_name="batch_1"

# build the data directory path using the variable
data_dir="./Corrected_Images/${batch_name}"
# get a list of all subdirectories in the Corrected_Images batch folder
mapfile -t image_dirs < <(find "$data_dir" -mindepth 1 -maxdepth 1 -type d | sed "s|^|$(pwd)/|")

echo "Number of image directories: ${#image_dirs[@]}"

for dir in "${image_dirs[@]}"; do
    echo "Found: $dir"
done

# loop over each directory and submit child jobs
for image_dir in "${image_dirs[@]}"; do
    # check job count for this user
    number_of_jobs=$(squeue -u "$USER" | wc -l)
    while [ "$number_of_jobs" -gt 990 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u "$USER" | wc -l)
    done
    sbatch cp_analysis_hpc_child.sh "$image_dir"
done

conda deactivate

echo "All CellProfiler jobs submitted!"
