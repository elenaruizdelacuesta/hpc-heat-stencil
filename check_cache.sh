#!/bin/bash
#SBATCH -A uTS25_Tornator_0
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH -t 00:05:00

lscpu | grep -i cache
lscpu 
