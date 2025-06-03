"""
This collection of functions runs CellProfiler in parallel and can convert the results into log files
for each process.
"""

import logging
import multiprocessing
import os
import pathlib
import subprocess
from concurrent.futures import Future, ProcessPoolExecutor
from typing import List

def run_cellprofiler_parallel(
    plate_info_dictionary: dict,
    run_name: str,
    group_level: str = "plate",
) -> None:
    """
    Run CellProfiler pipelines in parallel using multiprocessing.

    Args:
        plate_info_dictionary (dict): Dictionary of info to run CellProfiler. 
            If group_level="plate", keys are plate names.
            If group_level="well", keys are well names and each value must include a "plate_name".
        run_name (str): A name for this CellProfiler run (e.g., "whole_image_features").
        group_level (str): Indicates processing level, either "plate" or "well". Default is "plate".

    Raises:
        FileNotFoundError: If required paths don't exist.
        ValueError: If group_level is invalid.
        KeyError: If plate_name is missing when group_level is "well".
    """
    # Check if that group level is set properly
    if group_level not in {"plate", "well"}:
        raise ValueError("group_level must be either 'plate' or 'well'")

    commands = []
    log_dir = pathlib.Path("./logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    # Confirm that the plate info dictionary contains the necessary keys
    for name, info in plate_info_dictionary.items():
        if group_level == "well":
            if "plate_name" not in info:
                raise KeyError(
                    f"Missing 'plate_name' for well '{name}'",
                    "Please include 'plate_name' in each entry when running at the well level."
                )
            plate_name = info["plate_name"]
        else:
            plate_name = name

        # Set paths for pipeline and output
        path_to_pipeline = info["path_to_pipeline"]
        path_to_output = info["path_to_output"]
        pathlib.Path(path_to_output).mkdir(exist_ok=True)

        # Confirm that the pipeline exists
        if not pathlib.Path(path_to_pipeline).resolve(strict=True):
            raise FileNotFoundError(
                f"The file '{pathlib.Path(path_to_pipeline).name}' does not exist"
            )

        # If running with a loaddata CSV, run the correct command
        if "path_to_loaddata" in info:
            path_to_loaddata = info["path_to_loaddata"]
            command = [
                "cellprofiler",
                "-c",
                "-r",
                "-p",
                path_to_pipeline,
                "-o",
                path_to_output,
                "--data-file",
                path_to_loaddata,
            ]
        else:
            path_to_images = info["path_to_images"]
            if not pathlib.Path(path_to_images).is_dir():
                raise FileNotFoundError(
                    f"Directory '{pathlib.Path(path_to_images).name}' does not exist or is not a directory"
                )
            command = [
                "cellprofiler",
                "-c",
                "-r",
                "-p",
                path_to_pipeline,
                "-o",
                path_to_output,
                "-i",
                path_to_images,
            ]
            if group_level == "well":
                command.extend(["-g", f"Well={name}"])

        # Append the command to the commands list
        commands.append((command, name, plate_name))

    # Check CPU count
    cpu_count = os.cpu_count()
    if cpu_count is None:
        raise RuntimeError("Unable to determine CPU count with os.cpu_count(). Please check your system.")
    # Set maximum workers to CPU count minus one
    executor = ProcessPoolExecutor(max_workers=cpu_count - 1)
    futures = [
        executor.submit(subprocess.run, args=cmd, capture_output=True)
        for cmd, _, _ in commands
    ]

    # Collect results
    results = [future.result() for future in futures]
    print("All processes have been completed!")

    # Write results to log files
    for result, (_, name, plate_name) in zip(results, commands):
        if group_level == "plate":
            log_path = log_dir / f"{name}_{run_name}_run.log"
        else:
            sublog_dir = log_dir / plate_name
            sublog_dir.mkdir(parents=True, exist_ok=True)
            log_path = sublog_dir / f"{name}_{run_name}_run.log"

        with open(log_path, "w") as f:
            f.write(result.stderr.decode("utf-8"))

        if result.returncode != 0:
            print(
                f"A return code of {result.returncode} was returned for {name}, which means there was an error."
            )

    print("All results have been converted to log files!")

