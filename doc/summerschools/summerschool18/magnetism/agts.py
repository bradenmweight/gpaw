from myqueue.workflow import run


def workflow():
    return [
        task('xyz.py'),
        task('magnetism1.py', tmax='1h', deps='xyz.py'),
        task('magnetism2.py', tmax='2h', deps='xyz.py')]
