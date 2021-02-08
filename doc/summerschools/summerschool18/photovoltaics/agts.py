from myqueue.workflow import run


def workflow():
    return [
        task('pv1.py', tmax='13h'),
        task('pv2.py', tmax='13h', cores=8, deps='pv1.py'),
        task('pv3.py', tmax='13h', deps='pv2.py')]
