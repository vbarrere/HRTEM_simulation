#!/usr/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

np="${1:-1}"
export xyz_dir="/home/victor/Data/MD_Data/AgCo/Dataset1/XYZ"
export data_file="/home/victor/Data/MD_Data/AgCo/Dataset1/data.dat"

export substrate_dir="/home/victor/Data/carbon_substrates"
export output_data="dataset.npy" # descripteurs et images, un enregistrement par image (voir README.md)
export max_files=1 # Put -1 to use all files in the dataset
export atom_typ1="Ag"
export atom_typ2="Co"
export n_px=256
export nz=20
export ht=200.0
export n_rotations=1
export defocus_study=1 # 1 to study defocus and spherical aberration variation (Put 1 rotation per simulation)
export noise_study=0 # 1 to study noise variation (if 1, defocus_study must be 0 and if defocus_study is 1, noise_study must be 0)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 # un thread BLAS par rang MPI



start_time=$(date +%s)

mpirun --use-hwthread-cpus -np "$np" ./main

end_time=$(date +%s)
echo "Execution time: $((end_time - start_time)) seconds"
