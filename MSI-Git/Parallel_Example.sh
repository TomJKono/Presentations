#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=24,walltime=24:00:00
#PBS -A mcgaughs
#PBS -m abe
#PBS -M user@umn.edu
#PBS -q batch

# Load modules
module load parallel
module load clustalo

# Use this line for single-node parallelization
parallel < Parallel_Commands.txt

# Use these lines for multi-node parallelization. See
# https://www.msi.umn.edu/support/faq/how-can-i-use-gnu-parallel-run-lot-commands-parallel
# for more information. Note that the '--jobs 1' part is NOT USED, because
# each node is listed in the ${PBS_NODEFILE} file for each core that is
# allocated to the job.
export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY"
parallel --sshloginfile ${PBS_NODEFILE} --workdir ${PWD} < Parallel_Commands.txt
