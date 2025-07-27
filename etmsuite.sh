#!/bin/bash
#
#SBATCH --job-name=etmmoleculeRedox
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=23000MB 
#SBATCH --cpus-per-task=32
#SBATCH --time=120:00:00
#SBATCH --output=etmmoleculeRedox.out
#SBATCH --error=etmmoleculeRedox.err
#SBATCH -p qCPU120
#SBATCH -A CHEM9C4

source /sysapps/anaconda3/etc/profile.d/conda.sh
conda activate etmsuite

python main.py molecule.xyz -dns 0.5  -prop all -tpj 8 -ttt 32 -ki -exs 3


