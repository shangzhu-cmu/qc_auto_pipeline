from ase.db import connect
from gpaw import GPAW,Mixer,MixerDif,Davidson
import surf_autoconv as surf_ac
element='Li_mp-135'
element_bulk=connect('final_database/bulk.db').get(name=element)
h=element_bulk.h
xc=element_bulk.xc
sw=element_bulk.sw
spin=element_bulk.spin
#all settings but kpts (due to structure dependence)
struc='100'
calc=GPAW(xc=xc,
    h=h,
    symmetry = {'point_group': False},
    eigensolver=Davidson(3),
    mixer=Mixer(beta=0.05,nmaxold=5,weight=50),
    spinpol=spin,
    maxiter=333,
    occupations={'name': 'fermi-dirac','width': sw},
    poissonsolver={'dipolelayer': 'xy'})
surf_ac.surf_auto_conv(element,struc,calc,generator='pymatgen',init_layer=4,interval=2,fix_layer=2,vac=10,rela_tol=5,temp_print=True)