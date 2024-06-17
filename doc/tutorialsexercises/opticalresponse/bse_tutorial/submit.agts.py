from myqueue.workflow import run


def workflow():
    with run(script='gs_Si.py', cores=4, tmax='20m'):
        with run(script='eps_Si.py', cores=4, tmax='1h'):
            run(script='plot_Si.py')

    with run(script='gs_MoS2.py', cores=4, tmax='1h'):
        with run(script='pol_MoS2.py', cores=56, tmax='2h'):
            run(script='plot_MoS2.py')
        with run(script='get_2d_eps.py', tmax='2h'):
            run(script='plot_2d_eps.py')
            run(script='test_2d_eps.py')
        run(script='alpha_MoS2.py')
