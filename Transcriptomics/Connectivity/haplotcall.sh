#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=32 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=32 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=200G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=%J.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=%J.output ## stdout file name(The %J will print job ID number)
#SBATCH --mail-type=NONE ## Send your job status via e-mail: Valid type values are NONE, BEGIN, END, FAIL, REQUEUE, ALL
########### Job information #############
echo "================================"
echo "Start at `date`"
echo "Job id is $SLURM_JOBID"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NTASKS processors."
echo "================================"
#########################################

######## Load required modules ##########
#. /etc/profile.d/modules.sh # Required line for modules environment to work
#module load openmpi/1.8.4 python/2.7 # Load modules that are required by your program
#########################################
source /lustre1/home/mass/eskalon/miniconda/bin/activate gatk2
### Below you can enter your program job command ###
REF="GCA_052425735.1_crgOcupat_genomic.fasta"

TAGs=( \
op_A1 \
op_A3 \
op_B2 \
op_B5 \
op_C4 \
op_C6 \
op_D1 \
op_D3 \
op_E2 \
op_E5 \
op_F6 \
op_OO1 \
op_OO2 \
op_OO3 \
op_OO4 \
op_OO5 \
op_OO6 \
op_12 \
op_13 \
op_14 \
op_25 \
op_26 \
op_27 \
op_28 \
op_5 \
op_6 \
)
for TAG in ${TAGs[@]}; do

echo "${TAG} ${BASENAME} ${OUT}"

gatk --java-options "-Xmx200G -XX:ParallelGCThreads=32 -Djava.io.tmpdr=temp" HaplotypeCaller \
                --reference ${REF} \
        	--input ${TAG}.SplitNCigarReads.split.bam \
        	--output ${TAG}.HaplotypeCaller.g.vcf.gz \
        	-dont-use-soft-clipped-bases \
        	-ERC GVCF

done
