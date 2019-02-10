#!/bin/bash

# Utilisation exclusive du noeud
#MSUB -x
# Nb de noeuds
#MSUB -N 1
# Nb cores
#MSUB -c 32
# Temps max (s)
#MSUB -T 18000
# Outputs
#MSUB -o output/hasw_out.txt
#MSUB -e output/hasw_err.txt
# Partition
#MSUB -q haswell

# MSUB -r hybrid_Job

set -x
cd ${BRIDGE_MSUB_PWD}

### Load 'G++5'
module loadFile c++/gnu/5.2.0

# Set max thread number
export OMP_NUM_THREADS=32
export GOMP_CPU_AFFINITY=1-32

# Remove previous results
# Move to directory

# Executable
binary=./bin/bench_graph
nthreads=32
niter=1
bench=RMAT_ER
nb_v=16000000
nb_e=128000000


# Retrieve NUMA nodes stats
numastat > numa_stat.log

# Benchmark parameters
for ((i=1; $i <= $nthreads; i++)); do
  echo "perf record -g -s -F 99 $binary $bench $nthreads $niter $nb_v $nb_e" > out.log
  perf record -g -s -F 99 $binary $bench $i $niter $nb_v $nb_e
done
echo "Finished"
