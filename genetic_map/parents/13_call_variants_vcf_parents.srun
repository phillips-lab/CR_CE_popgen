#!/bin/bash
#SBATCH -A phillipslab
#SBATCH --partition=phillips       ### Partition
#SBATCH --job-name=parents    ### Job Name
#SBATCH --mem=20G
#SBATCH --time=50:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --array=0-5           ### Array index
#SBATCH --cpus-per-task=5            ### Number of CPU cores per task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#########################################################################
############# Variant calling for PX506 and PX553 ######################
#########################################################################
###!!! the previous steps of the analysis are scripts #1 - 8 in the diversity_stats/C. remanei



module load java easybuild GATK
ref="/projects/phillipslab/shared/ref_genomes/CR_PB_HIC/NCBI/CR.ncbi.softmasked.fasta"
tmp="/projects/phillipslab/ateterina/tmp"
GATK="/projects/phillipslab/ateterina/scripts/gatk-4.1.4.1/gatk"
cd /projects/phillipslab/ateterina/CR_popgen/data/reads/BAMS

declare -a CHROMS=("I" "II" "III" "IV" "V" "X")
CHR=${CHROMS[$SLURM_ARRAY_TASK_ID]}
echo $CHR
cd /projects/phillipslab/ateterina/CR_popgen/data/reads/BAMS

samples=$(find . | sed 's/.\///' | grep -E 'P(.)*.raw.new.g.vcf$' | sed 's/^/-V /')

#consolidate variants from many samples

$GATK --java-options "-Xmx15g -Xms10g" GenomicsDBImport -R $ref \
--genomicsdb-workspace-path Parents_database_${CHR} \
--tmp-dir $tmp -L ${CHR} --reader-threads 5 \
 $(echo $samples);




$GATK --java-options "-Xmx15g -Xms10g" GenotypeGVCFs \
       -R $ref \
       -V gendb://Parents_database_${CHR} \
       -O Parents_full_output.${CHR}.vcf \
       --include-non-variant-sites \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        -L ${CHR}
