if [ -n $PBS_ENVIRONMENT ]; then cd $PBS_O_WORKDIR; fi
GPAW=~/gpaw
ASE=~/ase
source /home/opt/modulefiles/modulefiles_el6.sh
module load GPAW
module load NUMPY/1.7.1-1
module load SCIPY/0.12.0-1
module load fftw
PLATFORM=linux-x86_64-$FYS_PLATFORM-2.6
export PATH=$GPAW/tools:$GPAW/build/bin.$PLATFORM:$PATH
export PYTHONPATH=$GPAW:$GPAW/build/lib.$PLATFORM:$PYTHONPATH
export PATH=$ASE/tools:$PATH
export PYTHONPATH=$ASE:$PYTHONPATH
alias gpaw-qsub=$GPAW/doc/platforms/Linux/Niflheim/gpaw-qsub.py