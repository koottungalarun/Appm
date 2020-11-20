#!/bin/bash

###################################################################################################
# Queue system requests

#SBATCH --job-name=appm		# job name displayed by squeue
###SBATCH --time=1:00:00 				# walltime, [HHH:MM:SS], abbreviated by -t
#SBATCH --nodes=1      				# number of cluster nodes, abbreviated by -N

#SBATCH -o slurm-%j.out 			# name of the stdout, using the job number (%j) and the first node (%N)

# for core spec. use either a fix number of cores with "ntasks" 
# or use just the pre defined 20 cores per node with "ntasks-per-node"
# with uncomenting "ntasks" you have just to spec. the number of nodes
###SBATCH --ntasks=40    			# number of MPI tasks, abbreviated by -n
#SBATCH --ntasks-per-node 1 		# each node has 20 cores

# additional information for allocated clusters
#SBATCH --account=915-0210     			# account - abbreviated by -A
###SBATCH --partition=  			# partition (queue), abbreviated by -p

###################################################################################################
# Basic setup

module purge
module add intelmpi/2018.1.163_64

HDF_INSTALL=/home/HSR/rfuchs/software/hdf5-1.12.0-linux-centos7-x86_64-shared

# ls 
echo $HDF_INSTALL
echo $MKLROOT

export LD_LIBRARY_PATH=$HDF_INSTALL/lib:$MKLROOT/lib/intel64;

echo "START"
./appm
wait 
echo "ENDING"
echo



