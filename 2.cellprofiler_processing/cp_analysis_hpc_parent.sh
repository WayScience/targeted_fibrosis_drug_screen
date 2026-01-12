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

# base directory for the corrected images
data_dir="./Corrected_Images/${batch_name}"

# count platemap folders (one level below batch folder)
mapfile -t platemap_dirs < <(find "$data_dir" -mindepth 1 -maxdepth 1 -type d)
echo "Number of platemap layouts in ${batch_name}: ${#platemap_dirs[@]}"

# get a list of all plates (folders starting with CARD) nested inside platemap layouts
mapfile -t plate_dirs < <(find "$data_dir" -mindepth 2 -maxdepth 2 -type d -name "CARD*" | sed "s|^|$(pwd)/|")
echo "Number of plates found: ${#plate_dirs[@]}"

# list plates found
for dir in "${plate_dirs[@]}"; do
    echo "Found plate: $dir"
done

# loop over each plate and submit child jobs
for plate_dir in "${plate_dirs[@]}"; do
    # check job count for this user
    number_of_jobs=$(squeue -u "$USER" | wc -l)
    while [ "$number_of_jobs" -gt 990 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u "$USER" | wc -l)
    done
    # submit the child script with the full plate path
    sbatch cp_analysis_hpc_child.sh --image_dir "$plate_dir"
done

conda deactivate

echo "All CellProfiler jobs submitted!"
