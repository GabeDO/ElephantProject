#! /bin/bash

#SBATCH --partition=Pegasus
#SBATCH --job-name=Slide
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8 
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=50:00:00

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""


cd $SLURM_SUBMIT_DIR
module load anaconda3
python SlidingWindow_DiversityMetrics.py "/nobackup/rogers_research/GabeO/Allignments/ZooElephants/SNPs/Elephant_SNPs_NoMissing.vcf.gz"


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

