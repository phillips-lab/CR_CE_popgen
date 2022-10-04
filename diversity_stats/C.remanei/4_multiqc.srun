#!/bin/bash
#SBATCH -A phillipslab
#SBATCH --partition=phillips       ### Partition
#SBATCH --job-name=mqc    ### Job Name
#SBATCH --time=00:240:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --array=0           ### Array index

module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 MultiQC/1.3-Python-3.6.1

cd projects/phillipslab/ateterina/CR_popgen/data/reads/fastqc_filt
multiqc .
