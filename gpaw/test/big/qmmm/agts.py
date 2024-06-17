from myqueue.workflow import run


def workflow():
    run(script='qmmm.py', cores=8, tmax='10m')
