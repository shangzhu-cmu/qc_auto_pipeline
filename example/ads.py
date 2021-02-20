import ads_selector as ads_sl
from gpaw import GPAW, MixerSum, Mixer, MixerDiff, Davidson
from ase.db import connect
element='Li_mp-46'
struc='100'
element_surf=connect('final_database/bulk.db').get(name=element+'('+struc+')')
h=element_surf.h
xc=element_surf.xc
sw=element_surf.sw
spin=element_surf.spin
kpts=[int(i) for i in (element_surf.kpts).split(',')]
calc=GPAW(xc=xc,
    h=h,
    kpts=kpts,
    symmetry = {'point_group': False},
    eigensolver=Davidson(3),
    mixer=Mixer(beta=0.05,nmaxold=5,weight=50),
    spinpol=spin,
    maxiter=333,
    occupations={'name': 'fermi-dirac','width': sw},
    poissonsolver={'dipolelayer': 'xy'})
ads_sl.ads_auto_select(element,
                struc,
                calc,
                ads='Li',
		        ads_pot_e=-1.89678,
                temp_print=True,
                size='1x1')