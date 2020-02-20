# flake8: noqa
from gpaw.test import equal
from ase import Atoms
from ase.units import Pascal, m
from ase.data.vdw import vdw_radii
from gpaw.mpi import rank
from gpaw import Mixer
from gpaw.solvation import (

def test_solvation_forces():
        SolvationGPAW,
        EffectivePotentialCavity,
        Power12Potential,
        LinearDielectric,
        KB51Volume,
        GradientSurface,
        VolumeInteraction,
        SurfaceInteraction,
        LeakedDensityInteraction)

    import numpy as np

    SKIP_ENERGY_CALCULATION = True
    F_max_err = 0.005

    h = 0.2
    u0 = 0.180
    epsinf = 80.
    T = 298.15
    atomic_radii = lambda atoms: [vdw_radii[n] for n in atoms.numbers]

    atoms = Atoms('NaCl', positions=((5.6, 5.6, 6.8), (5.6, 5.6, 8.8)))
    atoms.set_cell((11.2, 11.2, 14.4))


    atoms.calc = SolvationGPAW(
        mixer=Mixer(0.5, 7, 50.0),
        xc='oldPBE', h=h, setups={'Na': '1'},
        cavity=EffectivePotentialCavity(
            effective_potential=Power12Potential(atomic_radii, u0),
            temperature=T,
            volume_calculator=KB51Volume(),
            surface_calculator=GradientSurface()),
        dielectric=LinearDielectric(epsinf=epsinf),
        # parameters chosen to give ~ 1eV for each interaction
        interactions=[
            VolumeInteraction(pressure=-1e9 * Pascal),
            SurfaceInteraction(surface_tension=100. * 1e-3 * Pascal * m),
            LeakedDensityInteraction(voltage=10.)
        ]
    )


    def vac(atoms):
        return min(
            atoms.positions[0][2],
            14.4 - atoms.positions[1][2]
        )

    step = .05
    if not SKIP_ENERGY_CALCULATION:
        d = []
        E = []
        F = []
        while vac(atoms) >= 5.6:
            d.append(abs(atoms.positions[0][2] - atoms.positions[1][2]))
            E.append(atoms.calc.get_potential_energy(atoms, force_consistent=True))
            F.append(atoms.calc.get_forces(atoms))
            atoms.positions[0][2] -= step

        d = np.array(d)
        E = np.array(E)
        F = np.array(F)

        if rank == 0:
            np.save('d.npy', d)
            np.save('E.npy', E)
            np.save('F.npy', F)
            from pprint import pprint
            print('d')
            pprint(list(d))
            print()
            print('E')
            pprint(list(E))
            print()
            print('F')
            pprint([list([list(l2) for l2 in l1]) for l1 in F])
    else:
        # h=0.2, setups: 0.9.11271, analytic gradient for dielectric
        d = [
            2.0000000000000009,
            2.0500000000000007,
            2.1000000000000005,
            2.1500000000000004,
            2.2000000000000002,
            2.25,
            2.2999999999999998,
            2.3499999999999996,
            2.3999999999999995,
            2.4499999999999993,
            2.4999999999999991,
            2.5499999999999989,
            2.5999999999999988,
            2.6499999999999986,
            2.6999999999999984,
            2.7499999999999982,
            2.799999999999998,
            2.8499999999999979,
            2.8999999999999977,
            2.9499999999999975,
            2.9999999999999973,
            3.0499999999999972,
            3.099999999999997,
            3.1499999999999968,
            3.1999999999999966
        ]

        E = [
            -3.563293139400153,
            -3.8483728941723685,
            -4.0711353674406672,
            -4.2422702451403147,
            -4.370843143863639,
            -4.4644180141621437,
            -4.5291689227132208,
            -4.5703343350581118,
            -4.5923649740023134,
            -4.5988279473461287,
            -4.5927583650683008,
            -4.5766206693403326,
            -4.552512679729916,
            -4.5222435353897978,
            -4.4871564101028012,
            -4.4484418158342498,
            -4.4070815207527305,
            -4.3639583908891826,
            -4.3196547998536206,
            -4.2746896647162842,
            -4.2295530542966002,
            -4.1846071284157169,
            -4.1401332574837273,
            -4.0962878505444289,
            -4.0533187524468151
        ]

        F = [
            [[2.6695821634315027e-14, -1.9815222433664118e-14, -6.4048490688641513],
             [-3.4082241002543622e-12, -1.6240575950714694e-12, 6.4041543958163452]],
            [[5.0237752987734333e-14, -1.8243293058344225e-14, -5.0393351704101388],
             [-1.5274896963914418e-12, -2.6398128880388972e-12, 5.0384209277642569]],
            [[-3.3648148693779676e-14, -5.9189391352451007e-15, -3.9060170441670885],
             [-8.4063321767924475e-13, -1.8821508970704685e-12, 3.905543132953595]],
            [[8.033887557888106e-14, 2.2292480376612319e-14, -2.9689893614346223],
             [-1.223736631764757e-12, -3.0969447122690189e-12, 2.9698074624482889]],
            [[-7.1961182128217342e-14, -9.2369061233814621e-15, -2.1978946263288037],
             [-8.8389572960386298e-13, -1.7994648800211978e-12, 2.1973890330684451]],
            [[6.839819262929247e-14, -5.8813682508237254e-14, -1.5641714869085164],
             [-1.7622743497617748e-12, -3.0046379594525875e-12, 1.5635650982049627]],
            [[-9.2964056597828368e-15, -7.5214119452955559e-15, -1.0432778856648208],
             [-1.5739852469080523e-12, -2.4988123414408888e-12, 1.0427858515240367]],
            [[4.9365637419280049e-14, -1.2806637110762107e-14, -0.61777588826245566],
             [-4.4633495248199538e-12, -8.4857800635478483e-13, 0.6191354280714596]],
            [[-2.7465001491431585e-15, 2.2194992830262815e-15, -0.27335940186763846],
             [2.2159953452415821e-12, -1.6412130688619773e-12, 0.27444763375180131]],
            [[-1.1307851878181378e-15, 1.1152866935660596e-14, 0.0049642333143336626],
             [-4.8799705828538529e-13, -2.6066991522648012e-12, -0.0056198136590138066]],
            [[-3.347924840460358e-14, -3.2847274823028079e-14, 0.22997686507849791],
             [-1.3848503589740797e-12, -3.1093262854030672e-12, -0.23005715760038931]],
            [[-1.8200707826858258e-14, -4.3512573943996744e-14, 0.40941568333028211],
             [1.1291695481871819e-12, -1.7474358989696015e-12, -0.40852810117309907]],
            [[-1.8441970643900321e-14, 5.61993561755351e-15, 0.54932435973885352],
             [-2.9861492319230874e-12, -3.0227468023234651e-12, -0.54910391544690951]],
            [[5.7679940133610544e-14, -2.8084539456163582e-14, 0.65784317706542283],
             [-1.0234224221710884e-12, -2.5375102240498517e-12, -0.65870527795968759]],
            [[-4.036664664589713e-15, -4.8205738644069875e-15, 0.7418625322893313],
             [-1.2904684108015449e-12, -2.1348233844072512e-12, -0.74210985960596609]],
            [[1.5318647653685363e-14, -3.421407576047577e-14, 0.80443218412113771],
             [8.3926234818054389e-13, -2.6560910148144288e-12, -0.80377843042267316]],
            [[6.9157565324685627e-15, -2.0265394116808891e-14, 0.84757769602691035],
             [2.2564924910397676e-12, -2.2656092749597723e-12, -0.84788230013436461]],
            [[3.6816279559139879e-14, -5.6217408874846806e-14, 0.87630143837759911],
             [-2.8919485514440065e-12, -1.5785323374781975e-12, -0.8773089495263976]],
            [[-2.7265578221877411e-15, 1.8112097101058187e-14, 0.89457926961794421],
             [-4.3811750564077208e-12, -2.461714836405928e-12, -0.89468520749309732]],
            [[3.7482219205057949e-14, -4.1131307258530919e-14, 0.90318154899617831],
             [-1.0179375774172992e-12, -1.3810014419535958e-12, -0.90288363369205138]],
            [[6.8901452284609408e-14, -1.8368129697299518e-14, 0.90212377386324472],
             [-2.5725988597478254e-12, -2.4508096630572398e-12, -0.90278569263072728]],
            [[-2.9987853255786194e-14, -8.7933798499792121e-15, 0.89490384223117858],
             [-1.9656418164895607e-12, -2.497775366538778e-12, -0.89593617789571922]],
            [[3.5493492557854517e-14, -9.8132290416331934e-15, 0.88412548228462806],
             [-2.1672889098350696e-13, -1.9186968970762135e-12, -0.88434959210747588]],
            [[-4.2258218665175417e-14, -7.2197458334904237e-15, 0.86940669722572905],
             [-6.3029668473160416e-13, -2.8462694786350656e-12, -0.86883783222557032]],
            [[-1.3554138996860724e-14, 7.8037843761066284e-16, 0.84970515380081257],
             [-4.5784249812786949e-12, -2.0773781762455308e-12, -0.8501771195912331]]
        ]
        d = np.array(d)
        E = np.array(E)
        F = np.array(F)


    # test for orthogonal forces equal zero:
    equal(F[..., :2], .0, 1e-7)

    stencil = 2  # 1 is too rough, 3 does not change compared to 2
    FNa, FCl = F[..., 2].T
    FNa *= -1.
    # test symmetry
    equal(FNa, FCl, F_max_err)
    dd = np.diff(d)[0]
    kernel = {
        1: np.array((0.5, 0, -0.5)),
        2: np.array((-1. / 12., 2. / 3., 0, -2. / 3., 1. / 12.)),
        3: np.array((1. / 60., -0.15, 0.75, 0, -0.75, 0.15, -1. / 60.)),
    }

    dEdz = np.convolve(E, kernel[stencil] / step, 'valid')

    err = np.maximum(
        np.abs(-dEdz - FNa[stencil:-stencil]),
        np.abs(-dEdz - FCl[stencil:-stencil])
    )

    # test forces against -dE / dd finite difference
    print(err)
    equal(err, .0, F_max_err)

    if SKIP_ENERGY_CALCULATION:
        # check only selected points:

        def check(index):
            atoms.positions[0][2] = 6.8 - index * step
            F_check = atoms.get_forces()
            equal(F_check[..., :2], .0, 1e-7)
            FNa_check, FCl_check = F_check[..., 2].T
            FNa_check *= -1.
            equal(FNa_check, FCl_check, F_max_err)
            err = np.maximum(
                np.abs(-dEdz[index - stencil] - FNa_check),
                np.abs(-dEdz[index - stencil] - FCl_check)
            )
            print(err)
            equal(err, .0, F_max_err)

        l = len(FNa)
        # check(stencil)
        check(l // 2)
        # check(l - 1 - stencil)
