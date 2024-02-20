#!/bin/bash
#SBATCH -p dev_q
#SBATCH -A PALMORE_LG01
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH -t 2:00:00
#SBATCH -J test
#SBATCH -e %x.err
#SBATCH -o %x.out
#SBATCH --mail-user=linysh1997@vt.edu

# Load modules
source /home/linysh1997/.profile

mpirun -np $SLURM_NTASKS ./nga.dp.gnu.opt.mpi.exe -i input 

# We are Done
exit
