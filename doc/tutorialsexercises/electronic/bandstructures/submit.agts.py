from myqueue.workflow import run


def workflow():
    run(script='bandstructure.py')
    run(script='soc.py')
