""" TODO:

1. we should find a good way in which to store files elsewhere than static

2. currently the files that are not generated by weekly tests are copied
   from srcpath. This needs to be documented.

Make sure that downloaded files are copied to build dir on build
This must (probably) be done *after* compilation because otherwise dirs
may not exist.

"""
try:
    from urllib2 import urlopen, HTTPError
except ImportError:
    from urllib.request import urlopen
    from urllib.error import HTTPError
    import ssl
import os


srcpath = 'http://wiki.fysik.dtu.dk/gpaw-files'
agtspath = 'http://wiki.fysik.dtu.dk'


def get(path, names, target=None, source=None):
    """Get files from web-server.

    Returns True if something new was fetched."""

    if target is None:
        target = path
    if source is None:
        source = srcpath
    got_something = False
    # We get images etc from a web server with a self-signed certificate
    # That cause trouble on some machines.
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    for name in names:
        src = os.path.join(source, path, name)
        dst = os.path.join(target, name)

        if not os.path.isfile(dst):
            print(dst, end=' ')
            try:
                data = urlopen(src, context=ctx).read()
                with open(dst, 'wb') as sink:
                    sink.write(data)
                print('OK')
                got_something = True
            except HTTPError:
                print('HTTP Error!')
    return got_something


literature = """
askhl_10302_report.pdf  mortensen_gpaw-dev.pdf      rostgaard_master.pdf
askhl_master.pdf        mortensen_mini2003talk.pdf
marco_master.pdf        mortensen_paw.pdf           ss14.pdf
""".split()
get('doc/literature', literature, 'documentation')

# Note: bz-all.png is used both in an exercise and a tutorial.  Therefore
# we put it in the common dir so far, rather than any of the two places
get('.', ['bz-all.png'], 'static')

# These files have different destinations after webpage refactor then they
# used to have. Might need to keep track of this.
get('exercises/wavefunctions', ['co_bonding.jpg'],
    target='tutorialsexercises/wavefunctions/wavefunctions')
get('exercises/lrtddft', ['spectrum.png'],
    target='tutorialsexercises/opticalresponse/lrtddft')
get('tutorials/wannier90', ['GaAs.png', 'Cu.png', 'Fe.png'],
    target='tutorialsexercises/wavefunctions/wannier90')
get('tutorials/xas', ['h2o_xas_3.png', 'h2o_xas_4.png'],
    target='tutorialsexercises/opticalresponse/xas')
get('tutorials/xas',
    ['xas_illustration.png'], target='documentation/xas')

# This files is not used anymore?
# get('tutorialsexercises/opticalresponse/xas', ['xas_h2o_convergence.png'])
# ----

get('documentation/xc', 'g2test_pbe0.png  g2test_pbe.png  results.png'.split())
get('performance', 'dacapoperf.png  goldwire.png  gridperf.png'.split(),
    'static')

get('bgp', ['bgp_mapping_intranode.png',
            'bgp_mapping1.png',
            'bgp_mapping2.png'], 'platforms/BGP')

# workshop 2013 and 2016 photos:
get('workshop13', ['workshop13_01_33-1.jpg'], 'static')
get('workshop16', ['gpaw2016-photo.jpg'], 'static')


# files from http://wiki.fysik.dtu.dk/gpaw-files/things/

# Warning: for the moment dcdft runs are not run (files are static)!
dcdft_pbe_aims_stuff = """
dcdft_aims.tight.01.16.db.csv
dcdft_aims.tight.01.16.db_raw.csv
dcdft_aims.tight.01.16.db_Delta.txt
""".split()

get('things', dcdft_pbe_aims_stuff, target='setups')

# Warning: for the moment dcdft runs are not run (files are static)!
dcdft_pbe_gpaw_pw_stuff = """
dcdft_pbe_gpaw_pw.csv
dcdft_pbe_gpaw_pw_raw.csv
dcdft_pbe_gpaw_pw_Delta.txt
""".split()

get('things', dcdft_pbe_gpaw_pw_stuff, target='setups')

# Warning: for the moment dcdft runs are not run (files are static)!
dcdft_pbe_jacapo_stuff = """
dcdft_pbe_jacapo.csv
dcdft_pbe_jacapo_raw.csv
dcdft_pbe_jacapo_Delta.txt
""".split()

get('things', dcdft_pbe_jacapo_stuff, target='setups')

# Warning: for the moment dcdft runs are not run (files are static)!
dcdft_pbe_abinit_fhi_stuff = """
dcdft_pbe_abinit_fhi.csv
dcdft_pbe_abinit_fhi_raw.csv
dcdft_pbe_abinit_fhi_Delta.txt
""".split()

get('things', dcdft_pbe_abinit_fhi_stuff, target='setups')

g2_1_stuff = """
pbe_gpaw_nrel_ea_vs.csv pbe_gpaw_nrel_ea_vs.png
pbe_gpaw_nrel_opt_ea_vs.csv pbe_gpaw_nrel_opt_distance_vs.csv
pbe_nwchem_def2_qzvppd_opt_ea_vs.csv pbe_nwchem_def2_qzvppd_opt_distance_vs.csv
""".split()

get('things', g2_1_stuff, target='setups')

get('things', ['datasets.json'], 'setups')

# Carlsberg foundation figure:
get('.', ['carlsberg.png'])

get('static', ['NOMAD_Logo_supported_by.png'])

# Summer school 2022
get('summerschool2018',
    ['CreateTunnelWin.png', 'JupyterRunningMac.png', 'JupyterRunningWin.png',
     'Logged_in_Mac.png', 'Logged_in_Win.png', 'Moba_ssh.png',
     'UseTunnelWin.png'],
    target='summerschools/summerschool22')
get('summerschool2018',
    ['organometal.master.db'],
    target='summerschools/summerschool22/machinelearning')
get('summerschool2018',
    ['C144Li18.png', 'C64.png', 'final.png', 'initial.png',
     'Li2.png', 'lifepo4_wo_li.traj', 'NEB_init.traj'],
    target='summerschools/summerschool22/batteries')
get('summerschool2022',
    ['Intro_projects_CAMD2022.pdf'],
    target='summerschools/summerschool22')

def setup(app):
    # Get png and csv files and other stuff from the AGTS scripts that run
    # every weekend:
    from gpaw.doctools.agts_crontab import find_created_files

    for path in find_created_files():
        # the files are saved by the weekly tests under agtspath/agts-files
        # now we are copying them back to their original run directories
        if path.is_file():
            continue
        print(path, 'copied from', agtspath)
        get('agts-files', [path.name], str(path.parent), source=agtspath)
