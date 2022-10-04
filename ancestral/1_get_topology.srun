#!/bin/bash
#SBATCH -A phillipslab
#SBATCH --partition=phillips       ### Partition
#SBATCH --job-name=topology    ### Job Name
#SBATCH --time=100:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --array=0           ### Array index
#SBATCH --cpus-per-task=8            ### Number of CPU cores per task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK





#cactusdir='/projects/phillipslab/ateterina/scripts/progressiveCactus/bin'
#seqfile="CR.seqfile.txt"
workdir="/projects/phillipslab/ateterina/ANCESTRAL/CR/genomes"
chrname="X"


#C.remanei ancestral state reconstruction


# 1 step - extract 1 Mb in X chromosome from PX506
module load bedtools
#region X   10000000    11000000
#region I   9000000 10000000

bedtools getfasta -fi /projects/phillipslab/ateterina/ANCESTRAL/CR/genomes/GCA_010183535.1_CRPX506_genomic.chromosomes.fna -bed 1mb_region.${chrname}.bed -fo PX506.${chrname}.1Mb.fasta

# 2 step - blast in to other assemblies
#change the number of threads
module load easybuild ifort/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 NCBI-Toolkit/18.0.0

makeblastdb -in PX506.${chrname}.1Mb.fasta -parse_seqids -dbtype nucl



for seq in GCA_000149515.1_ASM14951v1_genomic.fna GCA_001643735.2_ASM164373v2_genomic.fna GCA_002259235.1_CaeLat1.0_genomic.fna GCA_002259225.1_CaeRem1.0_genomic.fna;do

    blastn -db PX506.${chrname}.1Mb.fasta -num_threads 8 -outfmt 6  -query $workdir/$seq -out ${seq/.fna/.blast_1Mb_${chrname}.out};

done



#extract those scaffolds from all other assemblies

###get a list of scaffolds

for i in *blast*_${chrname}.*out*;do
    awk '{if ($4>500) print;}' $i |cut -f1 - | uniq -c | sed -E "s/^( )*[0-9]+//" - > ${i/out/list.txt};
done

sed -i  's/ //g' *list.txt

###extract scaffolds

for seq in GCA_000149515.1_ASM14951v1_genomic.fna GCA_001643735.2_ASM164373v2_genomic.fna GCA_002259235.1_CaeLat1.0_genomic.fna GCA_002259225.1_CaeRem1.0_genomic.fna;do

    perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ${seq/fna/blast_1Mb_${chrname}.list.txt} $workdir/$seq > ${seq/.fna/.SUBSET.${chrname}.1Mb.fasta}

done

#run progressiveMAUVE

PX506="PX506.${chrname}.1Mb.fasta"
PX356="GCA_001643735.2_ASM164373v2_genomic.SUBSET.${chrname}.1Mb.fasta"
PB4641="GCA_000149515.1_ASM14951v1_genomic.SUBSET.${chrname}.1Mb.fasta"
CL="GCA_002259235.1_CaeLat1.0_genomic.SUBSET.${chrname}.1Mb.fasta"
PX439="GCA_002259225.1_CaeRem1.0_genomic.SUBSET.${chrname}.1Mb.fasta"

progressiveMauve="/projects/phillipslab/ateterina/scripts/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve"
out="CR_${chrname}_1Mb.pmauve"

$progressiveMauve --output=$out.xmfa --output-guide-tree=$out.tree --backbone-output=$out.backbone $PX506 $PX356 $PB4641 $PX439 $CL


#use it's tree fro progressiveCactus
