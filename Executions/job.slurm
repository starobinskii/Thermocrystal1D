#!/bin/bash
#SBATCH -n 1000
#SBATCH -J thCr1D
#SBATCH -o thCr1D-%j.out
#SBATCH -e thCr1D-%j.err
#SBATCH -t 30:00:00

mkdir -p Runs

for i in {1..50}; do
    sleep 1s
    time mpiexec.hydra ./Executions/task --path=./Runs --dispersion=1 --standing-wave --running-wave --save-velocities --env=cli --id=${i}
done
