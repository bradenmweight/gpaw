# web-page: h.png
def workflow():
    from myqueue.workflow import run
    with run(script='h.py'):
        run(function=plot)


def plot():
    import numpy as np
    import pylab as plt
    from ase.io import read
    with open('h.py').read().replace('ae', 'paw') as code:
        exec(code)
    ae = np.array([h.get_potential_energy() for h in read('H.ae.txt@:')])
    paw = np.array([h.get_potential_energy() for h in read('H.paw.txt@:')])
    ecut = range(200, 901, 100)
    plt.figure(figsize=(6, 4))
    plt.plot(ecut, ae[:-1] - ae[-1], label='ae')
    plt.plot(ecut, paw[:-1] - paw[-1], label='paw')
    plt.legend(loc='best')
    plt.xlabel('ecut [eV]')
    plt.ylabel('E(ecut)-E(1000 eV)')
    plt.savefig('h.png')
