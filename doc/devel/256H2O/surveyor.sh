qsub -A Gpaw -n 32 -t 40 --mode smp --env OMP_NUM_THREADS=1:GPAW_SETUP_PATH=$GPAW_SETUP_PATH:PYTHONPATH="${HOME}/gpaw:$PYTHONPATH":LD_LIBRARY_PATH="/bgsys/drivers/ppcfloor/gnu-linux/powerpc-bgp-linux/lib:$LD_LIBRARY_PATH" ${HOME}/gpaw/build/bin.linux-ppc64-2.5/gpaw-python ../b256H2O.py --sl_diagonalize=4,4,64,4 --gpaw=usenewlfc=1
