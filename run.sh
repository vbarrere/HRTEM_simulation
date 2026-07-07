#!/usr/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

np="${1:-48}"
export xyz_dir="../xyz2"
export images_data="images_full.dat"
export descriptors_data="descriptors_full.dat"
export max_files=-1 # Put -1 to use all files in the dataset
export atom_typ1="Ag"
export atom_typ2="Co"

start_time=$(date +%s)
mpirun --use-hwthread-cpus -np "$np" ./main

end_time=$(date +%s)
echo "Execution time: $((end_time - start_time)) seconds"
