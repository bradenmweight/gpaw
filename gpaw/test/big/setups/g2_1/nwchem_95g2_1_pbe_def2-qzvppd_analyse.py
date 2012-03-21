import os

import csv

import numpy as np

from ase.tasks.main import run
from ase.data.g2_1 import molecule_names
from ase.data.g2_1 import atom_names

ref = {'energy': {'NH3': -1537.6810863381743, 'S2': -21662.53340789671, 'SiH2_s3B1d': -7903.221529834168, 'Li': -203.05091256761085, 'CH3OH': -3146.3744968368742, 'SiH4': -7938.1877955148329, 'Si2H6': -15844.573176324142, 'PH3': -9333.1568491848029, 'PH2': -9315.9649519639133, 'HF': -2732.0171537716778, 'O2': -4088.6834866636477, 'SiH3': -7920.7014427048816, 'NH': -1501.3486099927698, 'Be': -398.09853859103947, 'SH2': -10863.760566567114, 'ClO': -14561.250898577431, 'H2O2': -4121.7339906521811, 'NO': -3532.3948238081593, 'ClF': -15231.878967564597, 'LiH': -218.96536083205314, 'HCO': -3095.9960806275412, 'CH3': -1082.5833067797735, 'CH4': -1100.9135366544626, 'Cl2': -25035.798062935362, 'Na': -4412.8227996807791, 'HOCl': -14578.847511658245, 'SiH2_s1A1d': -7903.9282286753705, 'SiO': -9920.0967334857633, 'F2': -5426.8218059748879, 'P2': -18569.620031115937, 'Si2': -15744.351053615766, 'C': -1028.5471953225028, 'CH': -1045.7495309245903, 'CO': -3081.3408011837355, 'CN': -2520.5041805686687, 'F': -2712.3065658470719, 'H': -13.60405569717515, 'LiF': -2921.2988166090859, 'O': -2041.248378488817, 'N': -1483.9814851078447, 'Na2': -8826.405269003686, 'P': -9282.2144876503189, 'Si': -7870.4445140239859, 'SO2': -14923.334564601748, 'NaCl': -16933.393059056176, 'Li2': -406.96307391293232, 'NH2': -1519.2263140507632, 'CS': -11864.981788520725, 'C2H6': -2169.2902518328342, 'N2': -2978.4659934121232, 'C2H4': -2135.8768034879554, 'HCN': -2540.108843934981, 'C2H2': -2102.0058009245777, 'CH2_s3B1d': -1064.0352818588406, 'CH3Cl': -13602.894948589868, 'BeH': -414.06750710089602, 'CO2': -5128.7887977210648, 'CH3SH': -11932.13562246209, 'OH': -2059.5496230411036, 'Cl': -12516.517962621117, 'S': -10828.82665694808, 'N2H4': -3041.6346796144649, 'H2O': -2078.4900285699109, 'SO': -12876.123333160007, 'CH2_s1A1d': -1063.3718592042937, 'H2CO': -3113.4561593161825, 'HCl': -12534.651514079676}, 'ae': {'NH3': 12.887434138803883, 'S2': 4.880094000549434, 'SiH2_s3B1d': 5.5689044158325487, 'CH3OH': 22.162700236853652, 'SiH4': 13.327058702147951, 'Si2H6': 22.059814093121531, 'PH3': 10.130194442959692, 'PH2': 6.5423529192448768, 'HF': 6.1065322274307618, 'O2': 6.1867296860136776, 'SiH3': 9.4447615893714101, 'NH': 3.7630691877498066, 'SH2': 7.725798224684695, 'ClO': 3.4845574674964155, 'H2O2': 12.029122280197043, 'NO': 7.1649602114976005, 'ClF': 3.0544390964078048, 'LiH': 2.3103925672671437, 'HCO': 12.596451119046378, 'CH3': 13.223944365744956, 'CH4': 17.950118543258895, 'Cl2': 2.7621376931274426, 'HOCl': 7.4771148511354113, 'SiH2_s1A1d': 6.2756032570350726, 'SiO': 8.4038409729600971, 'F2': 2.2086742807441624, 'P2': 5.1910558152994781, 'Si2': 3.4620255677946261, 'CH': 3.5982799049122605, 'CO': 11.545227372415866, 'CN': 7.9755001383209674, 'LiF': 5.9413381944032153, 'Na2': 0.7596696421278466, 'SO2': 12.011150676033139, 'NaCl': 4.052296754278359, 'Li2': 0.86124877771061392, 'NH2': 8.0367175485680491, 'CS': 7.607936250142302, 'C2H6': 30.571527004777181, 'N2': 10.503023196433787, 'C2H4': 24.366190054248818, 'HCN': 13.976107807457993, 'C2H2': 17.703298885221557, 'CH2_s3B1d': 8.2799751419872791, 'CH3Cl': 17.017623554724196, 'BeH': 2.3649128126813821, 'CO2': 17.744845420927959, 'CH3SH': 20.345547402808734, 'OH': 4.6971888551115626, 'N2H4': 19.255486610074513, 'H2O': 10.033538686743668, 'SO': 6.0482977231094992, 'CH2_s1A1d': 7.6165524874404582, 'H2CO': 16.452474110512412, 'HCl': 4.5294957613841689}}

atoms, task = run('nwchem molecule g2-1 -t 95g2_1_pbe_def2-qzvppd --atomize -s')

prop2molecules = {'ae': molecule_names, 'energy': molecule_names + atom_names}
prop2index = {'ae': 4, 'energy': 0}
prop2prec = {'ae': 0.01, 'energy': 0.01}

calc = {}

for p in ['ae', 'energy']:
    calc[p] = {}
    for m in prop2molecules[p]:
        try:
            calc[p][m] = task.results[m][prop2index[p]]
        except KeyError:
            calc[p][m] = None
            print 'Missing: ' + m + ' for property ' + p
            pass

    skeys = ref[p].keys()
    skeys.sort()
    for k in skeys:
        assert calc[p][k] is not None, 'Missing: ' + k + ' for property ' + p
        diff = calc[p][k] - ref[p][k]
        assert abs(diff) < prop2prec[p], 'Error: ' + k + ' ' + str(diff) + ' for property ' + p

    outfilename = 'nwchem_95g2_1_pbe_def2-qzvppd_%s.csv' % p

    d = []
    for k in skeys:
        d.append([k, calc[p][k]])
    csvwriter = csv.writer(open(outfilename, 'wb'))
    for row in d:
        csvwriter.writerow(row)
