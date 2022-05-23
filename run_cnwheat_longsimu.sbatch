#!/bin/bash
#SBATCH --account=inra_ea
#SBATCH --partition=defq
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=maurin.woussen@inrae.fr
#SBATCH --mail-type=ALL

# variable local au job
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# charge la bonne version de singularity
module load singularity/3.6.3

# run le main
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
	singularity run --writable-tmpfs cnwheat_longsimu_16-5-22.simg -l 2 -o cnwheat_caribu-active_append -n 2500 -s $SLURM_ARRAY_TASK_ID -w append

elif [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
	singularity run --writable-tmpfs cnwheat_longsimu_16-5-22.simg -l 2 -o cnwheat_ratp-active_append -n 2500 -s $SLURM_ARRAY_TASK_ID -w append

elif [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then
	singularity run --writable-tmpfs cnwheat_longsimu_16-5-22.simg -o cnwheat_vegetative_stages_default -n 2500 -s $SLURM_ARRAY_TASK_ID

elif [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then
	singularity run --writable-tmpfs cnwheat_longsimu_16-5-22.simg -l 2 -o cnwheat_ratp-active_final -n 2500 -s 2 -w final

fi
