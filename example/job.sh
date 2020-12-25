#!/bin/bash
#SBATCH -J test # Job name
#SBATCH -n 16 # Number of total cores
#SBATCH -N 1 # Number of nodes
#SBATCH --time=02-00:00
#SBATCH -A gpu
#SBATCH -p gpu
#SBATCH --mem-per-cpu=4000 # Memory pool for all cores in MB
#SBATCH -e test.err #change the name of the err file 
#SBATCH -o test.out # File to which STDOUT will be written %j is the job #
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kianpu20@cmu.edu # Email to which notifications will be sent

echo "Job started on `hostname` at `date`"

mpirun -np 16 gpaw python bulk_test.py

echo " "
echo "Job Ended at `date`"
