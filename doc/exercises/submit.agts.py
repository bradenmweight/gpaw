from myqueue.job import Job


def workflow():
    return [
        Job('water/h2o.py'),
        Job('wavefunctions/CO.py@8x15s'),
        Job('aluminium/Al_fcc.py@2x15s'),
        Job('aluminium/Al_bcc.py@2x15s'),
        Job('aluminium/Al_fcc_vs_bcc.py@2x15s'),
        Job('aluminium/Al_fcc_modified.py@2x15s'),
        Job('diffusion/initial.py@2x15s'),
        Job('diffusion/solution.py@2x15s'),
        Job('diffusion/densitydiff.py', deps=['diffusion/solution.py']),
        Job('wannier/si.py@8x15s'),
        Job('wannier/wannier-si.py', deps=['wannier/si.py']),
        Job('wannier/benzene.py@8x15s'),
        Job('wannier/wannier-benzene.py', deps=['wannier/benzene.py']),
        Job('band_structure/ag.py'),
        Job('vibrations/h2o.py@8x15s'),
        Job('vibrations/H2O_vib.py@8x15s', deps=['vibrations/h2o.py']),
        Job('vibrations/H2O_vib_2.py@4x15s', deps=['vibrations/H2O_vib.py']),
        Job('iron/ferro.py@4x15s'),
        Job('iron/anti.py@4x15s'),
        Job('iron/non.py@2x15s'),
        Job('iron/PBE.py',
            deps=['iron/ferro.py', 'iron/anti.py', 'iron/non.py']),
        Job('eels/test.py', deps=['band_structure/ag.py']),
        Job('gw/test.py'),
        Job('rpa/si.pbe.py'),
        Job('rpa/si.pbe+exx.py@4x15s', deps=['rpa/si.pbe.py']),
        Job('rpa/si.rpa_init_pbe.py'),
        Job('rpa/si.rpa.py@4x15s', deps=['rpa/si.rpa_init_pbe.py']),
        Job('stress/con_pw.py'),
        Job('stress/stress.py'),
        Job('transport/pt_h2_tb_transport.py'),
        Job('transport/pt_h2_lcao_manual.py'),
        Job('transport/pt_h2_lcao_transport.py',
            deps=['transport/pt_h2_lcao_manual.py'])]
