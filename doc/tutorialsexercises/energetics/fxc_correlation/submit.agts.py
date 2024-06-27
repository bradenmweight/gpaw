from myqueue.workflow import run


d_pbe = 'diamond.ralda_01_pbe.py'
d_corr = 'diamond.ralda_02_rapbe_rpa.py'
d_extra = 'diamond.ralda_03_extrapolate.py'

dk_lda = 'diam_kern.ralda_01_lda.py'
dk_dens = 'diam_kern.ralda_02_ralda_dens.py'
dk_wave = 'diam_kern.ralda_03_ralda_wave.py'
dk_rpa = 'diam_kern.ralda_04_rpa.py'
dk_comp = 'diam_kern.ralda_05_compare.py'


def workflow():
    with run(script='H.ralda_01_lda.py', cores=2, tmax='1m'):
        run(script='H.ralda_02_rpa_at_lda.py', cores=16, tmax='20m')
        run(script='H.ralda_03_ralda.py', cores=16, tmax='4h')
    with run(script='H.ralda_04_pbe.py', cores=2, tmax='1m'):
        run(script='H.ralda_05_rpa_at_pbe.py', cores=16, tmax='20m')
        run(script='H.ralda_06_rapbe.py', cores=16, tmax='4h')
    with run(script='CO.ralda_01_pbe_exx.py', cores=40, tmax='20m'):
        run(script='CO.ralda_02_CO_rapbe.py', cores=40, tmax='20h')
        c = run(script='CO.ralda_03_C_rapbe.py', cores=40, tmax='20h')
        o = run(script='CO.ralda_04_O_rapbe.py', cores=40, tmax='20h')
        with run(script=d_pbe, cores=16, tmax='1h'):
            with run(script=d_corr, cores=40, tmax='6h'):
                run(script=d_extra, cores=1, tmax='1m', deps=[c])
                run(script='CO.ralda_05_extrapolate.py',
                    deps=[c, o], cores=1, tmax='1m')
    with run(script=dk_lda, cores=8, tmax='1m'):
        with run(script=dk_dens, cores=40, tmax='4h'):
            with run(script=dk_wave, cores=8, tmax='1m'):
                with run(script=dk_rpa, cores=8, tmax='1m'):
                    run(script=dk_comp, cores=1, tmax='1m')
