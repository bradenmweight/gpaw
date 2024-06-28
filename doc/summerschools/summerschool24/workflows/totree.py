import asr
from ase.build import bulk


def materials():
    elements = ['Al', 'Si', 'Ti', 'Fe', 'Ni', 'Cu', 'Ag', 'Sb', 'Au']
    return {symbol: bulk(symbol) for symbol in elements}


workflow = asr.totree(materials(), name='material')
