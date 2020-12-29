from gpaw import GPAW
from ase.constraints import FixAtoms
from ase import build
from ase.io import write
from ase.db import connect
import os
import optimizer as opt
from ase.parallel import parprint
import numpy as np
import re

def surf_auto_conv(element,a,struc,num_of_layers=3,vac=8,fix_layer=2,h=0.14,k=6,xc='PBE',sw=0.1,rela_tol=10*10**(-3)):
    dispatcher = {'fcc100':build.fcc100,'fcc110':build.fcc110,'fcc111':build.fcc111,
                  'bcc100':build.bcc100,'bcc110':build.bcc110,'bcc111':build.bcc111,
                  }
    db_surf=connect(element+"_"+'surf'+'/'+'surf_converge.db')
    diff_primary=100
    diff_second=100
    iters=0
    #first optimize the layers
    while diff_primary>rela_tol and diff_second>rela_tol and iters < 6:
        slab = dispatcher[struc](element, size=(1, 1, num_of_layers),a=a,vacuum=vac)
        fix_mask=slab.positions[:,2] <= np.unique(slab.positions[:,2])[fix_layer-1]
        slab.set_constraint(FixAtoms(mask=fix_mask))
        slab.set_pbc([1,1,0])
        calc=GPAW(xc=xc,
            h=h,
            symmetry = {'point_group': False},
            kpts=(k,k,1),
            occupations={'name': 'fermi-dirac','width': sw},
            poissonsolver={'dipolelayer': 'xy'})
        slab.set_calculator(calc)
        name=element+"_"+"surf"+'/'+'results_layers'
        opt.surf_relax(slab, name, fmax=0.01, maxstep=0.04, replay_traj=None)
        db_surf.write(slab)
        layer=re.findall(r'\d+',str(slab.symbols))
        ß
    #second optimize the vaccume layer