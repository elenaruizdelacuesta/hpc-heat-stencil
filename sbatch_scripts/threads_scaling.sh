#!/bin/bash

echo "Threads scaling: single node with multiple threads"

# Configuración del problema
NODES=1
N_TASKS_PER_NODE=1
TOTAL_TASKS=1
N_STEPS=750
GRID_SIZE_X=15000
GRID_SIZE_Y=15000

# Serie de hilos a probar
THREADS_LIST=(1 2 4 8 16 32 56 84 112)

for OMP_THREADS in "${THREADS_LIST[@]}"; do
    JOB_NAME="thread_scaling_${OMP_THREADS}_threads"
    echo "Submitting job ${JOB_NAME} with ${OMP_THREADS} threads"

    sbatch \
        --export=ALL,GRID_SIZE_X=${GRID_SIZE_X},GRID_SIZE_Y=${GRID_SIZE_Y},N_STEPS=${N_STEPS},OMP_THREADS=${OMP_THREADS},JOB_NAME=${JOB_NAME},TOTAL_TASKS=${TOTAL_TASKS} \
        --nodes=${NODES} \
        --ntasks-per-node=${N_TASKS_PER_NODE} \
        --cpus-per-task=${OMP_THREADS} \
        --job-name=${JOB_NAME} \
        sbatch_scripts/go_dcgp.sh
done

echo "All threads scaling jobs submitted"
