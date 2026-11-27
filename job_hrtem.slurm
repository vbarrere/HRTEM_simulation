#!/bin/bash

#SBATCH --job-name=hrtem_sim           # label in squeue; becomes %x in the log name
#SBATCH --account=rsk@cpu              # allocation charged; also routes to the V100 nodes, i.e. GPUs, change to "rsk@cpu" in case of CPU
#SBATCH --nodes=8                     # one machine GPU, in case of CPU 2 CPU with 20 cores, ie 40 threads, without hyperthreading
#SBATCH --ntasks=320                  # one process (a single python); in case of CPU set to the number of threads you want; e.g. 40 in case of 1 node, or 80 for two nodes, or 80 for one core and hyperthreading enabled
##SBATCH --gres=gpu:1                  # one GPU; a node has four ; remove this line for a CPU only job
#SBATCH --cpus-per-task=1              # 10 physical cores = 1/4 node = 40 GB RAM; keep 10 per GPU, set to 1 in case of cpu
#SBATCH --hint=nomultithread           # count physical cores, not hardware threads; otherwise "--hint=multithread" for the use of hyperthreading
#SBATCH --time=01:00:00                # hh:mm:ss; 10 min default; max 20 h (100 h with qos_gpu-t4)
#SBATCH --output=slurm_%x_%j.out       # stdout and stderr together; %j = job id
#SBATCH --qos=qos_cpu-dev
##SBATCH --qos=qos_gpu-dev             # tests up to 2 h, starts faster
##SBATCH --qos=qos_gpu-t4              # jobs longer than 20 h, up to 100 h
##SBATCH --constraint=v100-32g         # only if 16 GB of GPU memory is not enough

# --- environment ------------------------------------------------------------
cd "${SLURM_SUBMIT_DIR}"                       # where sbatch was run; relative paths below are from here
source /etc/profile.d/modules.sh   2>/dev/null || true
source /etc/profile.d/z_modules.sh 2>/dev/null || true
module purge                                   # ignore whatever the login shell had loaded
module load netcdf-fortran/4.5.3-mpi-cuda
module load fftw/3.3.10-mpi-cuda
module load openblas/0.3.10
#export OMP_NUM_THREADS=1               # one thread per rank, no oversubscription, for lammps

# --- bookkeeping in the log ---------------------------------------------------
echo "job ${SLURM_JOB_ID} on $(hostname), started $(date)"
#nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true  # write details of the GPU to log

# --- the work -----------------------------------------------------------------
# Everything read here must already be on Jean Zay: compute nodes have no internet.
export xyz_dir="/lustre/fsn1/projects/rech/rsk/usr18vl/Data/MD_Data/Dataset1/XYZ"
export data_file="/lustre/fsn1/projects/rech/rsk/usr18vl/Data/MD_Data/Dataset1/data.dat"
export substrate_dir="/lustre/fsn1/projects/rech/rsk/usr18vl/Data/carbon_substrates"
export output_data="data.npy"
export max_files=-1 # Put -1 to use all files in the dataset
export atom_typ1="Ag"
export atom_typ2="Co"
export n_px=256
export nz=20
export ht=200.0
export n_rotations=10 # Put 1 for defocus_study, 10 otherwise
export defocus_study=0 # 1 to study defocus and spherical aberration variation (Put 1 rotation per simulation)
export noise_study=0 # 1 to study noise variation (if 1, defocus_study must be 0 and if defocus_study is 1, noise_study must be 0)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 # un thread BLAS par rang MPI

make

start_time=$(date +%s)

srun ./main
rc=$?
#mpirun --use-hwthread-cpus -np "640" ./main

end_time=$(date +%s)
echo "Execution time: $((end_time - start_time)) seconds"

#srun lmp -in lammps_script             # never mpirun -np 40; srun already knows the 40 ranks

echo "finished $(date) with exit code ${rc}"

# ~10 min for one subdataset, ~40min-1h for the full dataset, on 16 nodes with 40 threads each (640 ranks) (use 8 nodes and 320 ranks for a faster job, but less parallelism ~15-20 min for one subdataset)
# ~45 min for defocus analysis (put only 1 rotation)
# ~5 min for noise analysis