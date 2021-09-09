#!/bin/bash
#SBATCH --job-name=AB_GIC_run12
#SBATCH --time=7-00:00
#SBATCH --mem-per-cpu=4000M
#SBATCH --ntasks=71
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dcordell@ualberta.ca
#SBATCH --account=rrg-unsworth-ab

cd $SLURM_SUBMIT_DIR
echo "prog started at: `date`"
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo ""
srun /home/dcordell/3D/f90/Mod3DMT -F LAB_sediments_topo_SABC_deep_NEcraton.model modem_topo.data fwd.dat EMsol control.fwd
echo "prog finished at: `date`"