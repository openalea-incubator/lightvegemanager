#!/bin/bash
#SBATCH --account=inra_ea
#SBATCH --partition=defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maurin.woussen@inrae.fr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1


# variable local au job
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# charge la bonne version de singularity
module load singularity/3.6.3

# run le main
singularity run analyse_tesselation_3-5-22.simg -l 0 -o /nfs/work/inra_ea/P3F/mwoussen/analyse_tesselation/ -i 64

# crée un zip des résultats
# cd /nfs/work/inra_ea/P3F/mwoussen
# zip -r analyse_tesselation_test.zip analyse_tesselation

# affiche heure de fin
date
