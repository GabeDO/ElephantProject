#! /bin/bash

#SBATCH --partition=Pegasus
#SBATCH --job-name=PCA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
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
echo "Im makin snakes on on compute node $(/bin/hostname -s)"


module load anaconda3
module load R/4.1.3

cd ..

echo "python PythonScripts/FISandHETp10KB.py"
Rscript --vanilla PythonScripts/MakeHetGraph.R


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"


