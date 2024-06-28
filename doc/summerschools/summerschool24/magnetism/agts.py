from myqueue.workflow import run


def workflow():
    with run(script='xyz.py'):
        run(script='magnetism1.py', tmax='1h')
        run(script='magnetism2.py', tmax='2h')
