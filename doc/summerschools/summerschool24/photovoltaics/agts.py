from myqueue.workflow import run


def workflow():
    with run(script='pv1.py', tmax='13h'):
        with run(script='pv2.py', tmax='13h', cores=8):
            run(script='pv3.py', tmax='13h')
