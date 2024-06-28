import os
import shutil
from pathlib import Path
from myqueue.workflow import run


def workflow():
    if os.getenv('AGTS_FILES'):
        dir = Path(os.getenv('AGTS_FILES'))
        file = Path('organometal.db')
        if not file.is_file():
            shutil.copyfile(dir / file, file)
    run(script='machinelearning.py', tmax='8h')
