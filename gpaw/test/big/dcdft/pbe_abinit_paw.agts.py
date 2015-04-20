# Note: due to how agts works external executables (abinit, aims, etc.)
# must be run on the submission platform (bmaster1), currently opteron4
def agts(queue):
    queue.add('pbe_abinit_paw.py Au', ncpus=1,
              queueopts='-l nodes=1:ppn=4:opteron4', walltime=10*60)