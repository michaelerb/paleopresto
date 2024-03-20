#!/bin/bash
#SBATCH --job-name=proxy_dashboards  # Name of the job
#SBATCH --output=output.txt          # File for output and errors
#SBATCH --time=1:00:00               # Maximum time for job to run
#SBATCH --mem=10000                  # Memory (MB)

# Run this with the command: sbatch run_script.sh.
srun python -u make_proxy_images_1_temp12k.py '12ka'
srun python -u make_proxy_images_1_temp12k.py '21ka'
srun python -u make_proxy_images_1_temp12k.py 'all'

