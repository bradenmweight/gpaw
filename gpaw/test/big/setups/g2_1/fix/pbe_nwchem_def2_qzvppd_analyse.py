import os

import csv

import numpy as np

from ase.tasks.main import run
from ase.data.g2_1 import molecule_names
from ase.data.g2_1 import atom_names

# these results are calculated at g2_1 geometries

ref = {'energy': {'NH3': -1537.8932102701845, 'S2': -21662.65237345071, 'SiH2_s3B1d': -7903.3651685467157, 'Li': -203.05091256761085, 'CH3OH': -3146.7614314384177, 'SiH4': -7938.4615430723688, 'Si2H6': -15845.077065956493, 'PH3': -9333.3945046039225, 'PH2': -9316.1237401236212, 'HF': -2732.0738778569244, 'O2': -4088.7018992526946, 'SiH3': -7920.9077574897519, 'NH': -1501.4239187420515, 'Be': -398.09853859103947, 'SH2': -10863.93376081882, 'ClO': -14561.285408235011, 'H2O2': -4121.9417690700211, 'NO': -3532.6816767356731, 'ClF': -15231.960222248297, 'LiH': -218.96907143369367, 'HCO': -3096.1951818382831, 'CH3': -1082.7988715085448, 'CH4': -1101.1756308328961, 'Cl2': -25035.885491471996, 'Na': -4412.8227996807791, 'HOCl': -14578.980544432243, 'SiH2_s1A1d': -7904.0658697482577, 'SiO': -9920.205514840487, 'F2': -5426.9000090032869, 'P2': -18569.696505717548, 'Si2': -15744.414966665674, 'C': -1028.5471958466758, 'CH': -1045.8215809606679, 'CO': -3081.4426304822496, 'CN': -2521.0068098554066, 'F': -2712.3065640405075, 'H': -13.60405569717515, 'LiF': -2921.366271359771, 'O': -2041.2483247580117, 'N': -1483.9814851078447, 'Na2': -8826.4081399044881, 'P': -9282.2144876503189, 'Si': -7870.4445140227344, 'SO2': -14923.470837520383, 'NaCl': -16933.432214428689, 'Li2': -406.96848818968698, 'NH2': -1519.3729803245799, 'CS': -11865.161101293274, 'C2H6': -2169.7958129976828, 'N2': -2978.476002513908, 'C2H4': -2136.2925745709285, 'HCN': -2540.2599155377502, 'C2H2': -2102.2847792279099, 'CH2_s3B1d': -1064.1843674510139, 'CH3Cl': -13603.219836063943, 'BeH': -414.10792412473955, 'CO2': -5129.0801085274989, 'CH3SH': -11932.531761054999, 'OH': -2059.6210065114869, 'Cl': -12516.517962621823, 'S': -10828.826656943618, 'N2H4': -3042.0310600740677, 'H2O': -2078.6212440372983, 'SO': -12876.189938295027, 'CH2_s1A1d': -1063.5112436284683, 'H2CO': -3113.7294275880872, 'HCl': -12534.741344990156}, 'ea': {'NH3': 13.099558070814055, 'S2': 4.9990595634735655, 'SiH2_s3B1d': 5.7125431296317402, 'CH3OH': 22.549688045029598, 'SiH4': 13.600806260935315, 'Si2H6': 22.563703727975735, 'PH3': 10.367849862079311, 'PH2': 6.7011410789527872, 'HF': 6.1632581192416183, 'O2': 6.2052497366712487, 'SiH3': 9.6510763754931759, 'NH': 3.8383779370315096, 'SH2': 7.89899248085203, 'ClO': 3.5191208551768796, 'H2O2': 12.237008159648212, 'NO': 7.4518668698165129, 'ClF': 3.1356955859664595, 'LiH': 2.3141031689076783, 'HCO': 12.79560553642068, 'CH3': 13.439508570343378, 'CH4': 18.21221219751942, 'Cl2': 2.8495662283494312, 'HOCl': 7.6102013552317658, 'SiH2_s1A1d': 6.4132443311737006, 'SiO': 8.5126760597413522, 'F2': 2.2868809222718482, 'P2': 5.2675304169097217, 'Si2': 3.5259386202051246, 'CH': 3.6703294168169123, 'CO': 11.647109877562343, 'CN': 8.4781289008860767, 'LiF': 6.0087947516526583, 'Na2': 0.76254054292985529, 'SO2': 12.14753106074204, 'NaCl': 4.0914521260856418, 'Li2': 0.86666305446527758, 'NH2': 8.1833838223847124, 'CS': 7.7872485029802192, 'C2H6': 31.077087121279874, 'N2': 10.513032298218604, 'C2H4': 24.781960088876076, 'HCN': 14.12717888605448, 'C2H2': 17.982276140207887, 'CH2_s3B1d': 8.4290602099877106, 'CH3Cl': 17.342510503920494, 'BeH': 2.4053298365249134, 'CO2': 18.036263164799493, 'CH3SH': 20.741685476006751, 'OH': 4.768626056300036, 'N2H4': 19.651867069677337, 'H2O': 10.164807884936181, 'SO': 6.1149565933974372, 'CH2_s1A1d': 7.7559363874420342, 'H2CO': 16.72579558904954, 'HCl': 4.6193266711579781}}

tag = 'pbe_nwchem_def2_qzvppd'

def main():

    atoms, task = run('nwchem molecule g2_1 -t ' + tag + ' --atomize -s')

    prop2molecules = {'ea': molecule_names, 'energy': molecule_names + atom_names}
    prop2index = {'ea': 'atomic energy', 'energy': 'energy'}
    prop2prec = {'ea': 0.01, 'energy': 0.01}

    calc = {}

    for p in ['energy', 'ea']:
        calc[p] = {}
        for m in prop2molecules[p]:
            try:
                if p == 'energy':
                    calc[p][m] = task.data[m][prop2index[p]]
                else:
                    # atomization energy
                    calc[p][m] = - calc['energy'][m]
                    calc[p][m] += task.data[m][prop2index[p]] # energy of atoms
            except (KeyError, TypeError):
                calc[p][m] = None
                print 'Missing: ' + m + ' for property ' + p
                pass

        skeys = ref[p].keys()
        skeys.sort()
        for k in skeys:
            assert calc[p][k] is not None, 'Missing: ' + k + ' for property ' + p
            diff = calc[p][k] - ref[p][k]
            assert abs(diff) < prop2prec[p], 'Error: ' + k + ' ' + str(diff) + ' for property ' + p

        outfilename = tag + '_%s.csv' % p

        d = []
        for k in skeys:
            d.append([k, calc[p][k]])
        csvwriter = csv.writer(open(outfilename, 'wb'))
        for row in d:
            csvwriter.writerow(row)

if __name__ == '__main__':
    main()
