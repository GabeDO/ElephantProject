#! /bin/bash

#SBATCH --partition=Pegasus
#SBATCH --job-name=FF_SnakeMake_Pipeline
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=200:00:00


echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

cd $SLURM_SUBMIT_DIR
echo "Im makin snakes on compute node $(/bin/hostname -s)"

module load samtools
module load bwa
module load snakemake
module load bedtools2
module load blast
module load plink

#snakemake -j 40 --use-conda --rerun-incomplete

snakemake --slurm  -p --profile Config/ --use-conda --rerun-incomplete --latency-wait 60 --jobs $SLURM_NTASKS --cores $SLURM_NTASKS


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
