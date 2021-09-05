#!/bin/bash
##############################################################
#                    SIMULATION sbatch file                  #
#                                                            #
#                                                            #
#               -Put vigilance stage NREM or WAKE,           #
#		          and the root directory as an               #
#                   argument to sbatch file                  #
##############################################################
#SBATCH --time=0-05:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5200M
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=farhad.razi01@estudiant.upf.edu
#SBATCH --output=output_run_%x.out
#SBATCH --error=error_run_%x.err
# 
# 
# 
# 
ml Python
# 
# 
# 
# 
export VIG_STATE=$1 ROOT=$2 ROOT_C=$3
# 
# 
# 
# 
mpirun -np $SLURM_NTASKS python SIMULATION.py
