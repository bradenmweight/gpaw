from myqueue.workflow import run


def workflow():
    return [
        task('calculate.py@1:1h'),
        task('plot.py', deps='calculate.py')]
