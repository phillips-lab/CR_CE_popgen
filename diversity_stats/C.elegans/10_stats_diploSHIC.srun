#!/bin/bash
#SBATCH -A phillipslab
#SBATCH --partition=phillips       ### Partition
#SBATCH --job-name=12statsCE    ### Job Name
#SBATCH --nodelist=n219
#SBATCH --time=15:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --array=0           ### Array index
#SBATCH --cpus-per-task=1            ### Number of CPU cores per task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


#module load java

module load java easybuild GATK bedtools samtools


diploSHIC="/projects/phillipslab/ateterina/scripts/diploSHIC/diploSHIC.py"
ref="/projects/phillipslab/ateterina/CE_haw_subset/ref_245/c_elegans.PRJNA13758.WS245.genomic.fa"
VCF="CE_population_filt_snps_5-100_0.5_fin.vcf"
MASK="CE_mask_these_region5-100_0.5.bed"
name=${MASK/_these_region/}



cd /projects/phillipslab/ateterina/CE_haw_subset/data/BAMS

bedtools maskfasta -fi $ref -bed $MASK -fo ${name/.bed/.fasta}
samtools faidx ${name/.bed/.fasta}

#I	15072434	3	50	51
#II	15279421	15373890	50	51
#III	13783801	30958905	50	51
#IV	17493829	45018387	50	51
#V	20924180	62862096	50	51
#X	17718942	84204763	50	51


#not normalized stats!
#a sligtly modified version of diploSHIC, see the "diploSHIC_note.txt" for detailes, it estimates non normaalized stats and 1 subwindow per window

python $diploSHIC fvecVcf --winSize 100000 --maskFileName ${name/.bed/.fasta} --unmaskedFracCutoff 0.10 --numSubWins 1 diploid $VCF I 15072434 ${VCF/.vcf/}.NOPHASE.I.100K.0.1.stats
python $diploSHIC fvecVcf --winSize 100000 --maskFileName ${name/.bed/.fasta} --unmaskedFracCutoff 0.10 --numSubWins 1 diploid $VCF II 15279421 ${VCF/.vcf/}.NOPHASE.II.100K.0.1.stats
python $diploSHIC fvecVcf --winSize 100000 --maskFileName ${name/.bed/.fasta} --unmaskedFracCutoff 0.10 --numSubWins 1 diploid $VCF III 13783801 ${VCF/.vcf/}.NOPHASE.III.100K.0.1.stats
python $diploSHIC fvecVcf --winSize 100000 --maskFileName ${name/.bed/.fasta} --unmaskedFracCutoff 0.10 --numSubWins 1 diploid $VCF IV 17493829 ${VCF/.vcf/}.NOPHASE.IV.100K.0.1.stats
python $diploSHIC fvecVcf --winSize 100000 --maskFileName ${name/.bed/.fasta} --unmaskedFracCutoff 0.10 --numSubWins 1 diploid $VCF V 20924180 ${VCF/.vcf/}.NOPHASE.V.100K.0.1.stats
python $diploSHIC fvecVcf --winSize 100000 --maskFileName ${name/.bed/.fasta} --unmaskedFracCutoff 0.10 --numSubWins 1 diploid $VCF X 17718942 ${VCF/.vcf/}.NOPHASE.X.100K.0.1.stats
