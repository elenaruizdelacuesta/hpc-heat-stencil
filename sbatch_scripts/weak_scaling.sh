#!/bin/bash

echo "Weak scaling: multinode study with constant workload per task"

N_STEPS=500
TASKS_PER_NODE=8
OMP_THREADS=14
CPUS_PER_TASK=${OMP_THREADS}
LOCAL_SIZE=4000   # tamaño de subgrilla por tarea

for NODES in 1 2 4 8 16; do
    TOTAL_TASKS=$((NODES * TASKS_PER_NODE))

    case ${TOTAL_TASKS} in
        8) px=4; py=2 ;;
       16) px=4; py=4 ;;
       32) px=8; py=4 ;;
       64) px=8; py=8 ;;
      128) px=16; py=8 ;;
    esac

    GRID_SIZE_X=$((px * LOCAL_SIZE))
    GRID_SIZE_Y=$((py * LOCAL_SIZE))
    
    JOB_NAME="weak_scale_${NODES}n_${TOTAL_TASKS}t"

    echo "Submitting ${JOB_NAME}: ${NODES} nodes, ${TOTAL_TASKS} tasks, grid ${GRID_SIZE_X}x${GRID_SIZE_Y}"

    sbatch \
        --nodes=${NODES} \
        --ntasks=${TOTAL_TASKS} \
        --ntasks-per-node=${TASKS_PER_NODE} \
        --cpus-per-task=${CPUS_PER_TASK} \
        --job-name=${JOB_NAME} \
        --export=ALL,GRID_SIZE_X=${GRID_SIZE_X},GRID_SIZE_Y=${GRID_SIZE_Y},N_STEPS=${N_STEPS},OMP_THREADS=${OMP_THREADS},TOTAL_TASKS=${TOTAL_TASKS},JOB_NAME=${JOB_NAME} \
        sbatch_scripts/go_dcgp.sh
done

echo "All weak scaling jobs submitted"
