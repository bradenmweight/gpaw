import runpy
globals().update(runpy.run_path('workflow.py'))


@asr.parametrize_glob('*/material')
def workflow(material):
    calculator = {'mode': 'pw', 'kpts': {'density': 1.0}, 'txt': 'gpaw.txt'}
    return MyWorkflow(atoms=material, calculator=calculator)
