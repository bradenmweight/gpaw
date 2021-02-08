from myqueue.workflow import run


def workflow():
    return [
        task('gs_3x3_defect.py@16:5m'),
        task('unfold_3x3_defect.py@16:10m', deps='gs_3x3_defect.py'),
        task('plot_sf.py', deps='unfold_3x3_defect.py')]
