from gpaw import GPAW
from ase.build import bulk
from ase.io import write
from ase.db import connect
import os
import optimizer as opt
from ase.parallel import parprint
import numpy as np

def bulk_auto_conv(element,h=0.16,k=6,xc='PBE',sw=0.1,rela_tol=10*10**(-3)):
    ## TO-DO: need a better way to make directory without outputting error
    #os.makedirs(element+"_"+'bulk')
    db_h=connect(element+"_"+'bulk'+'/'+'grid_converge.db')
    db_k=connect(element+"_"+'bulk'+'/'+'kpts_converge.db')
    db_sw=connect(element+"_"+'bulk'+'/'+'sw_converge.db')
    diff_primary=100
    diff_second=100
    iters=0
    #start with grid spacing convergence
    while diff_primary>rela_tol and diff_second>rela_tol and iters < 6:
        if iters>0:
            h=np.round(h-0.02,decimals=2)
        atoms=bulk(element)
        calc=GPAW(xc=xc,h=h,kpts=(k,k,k),occupations={'name':'fermi-dirac','width':sw})
        atoms.set_calculator(calc)
        opt.optimize_bulk(atoms,step=0.05,fmax=0.01,location=element+"_"+'bulk'+'/'+'results_grid',extname='{}'.format(h))
        db_h.write(atoms,h=h)
        if iters>=2:
            fst=db_h.get_atoms(id=iters-1)
            snd=db_h.get_atoms(id=iters)
            trd=db_h.get_atoms(id=iters+1)
            diff_primary=max(abs(snd.get_potential_energy()-fst.get_potential_energy()),
                            abs(trd.get_potential_energy()-fst.get_potential_energy()))
            diff_second=abs(trd.get_potential_energy()-snd.get_potential_energy())
        iters+=1

    #kpts convergence
    diff_primary=100
    diff_second=100
    iters=0
    while diff_primary>rela_tol and diff_second>rela_tol and iters<6: 
        if iters>0:
            k=int(k+2)
        atoms=bulk(element)
        calc=GPAW(xc=xc,h=h,kpts=(k,k,k),occupations={'name':'fermi-dirac','width':sw})
        atoms.set_calculator(calc)
        opt.optimize_bulk(atoms,step=0.05,fmax=0.01,location=element+"_"+'bulk'+'/'+'results_kpts',extname='{}'.format(k))
        db_k.write(atoms,kpts=k)
        if iters>=2:
            fst=db_h.get_atoms(id=iters-1)
            snd=db_h.get_atoms(id=iters)
            trd=db_h.get_atoms(id=iters+1)
            diff_primary=max(abs(snd.get_potential_energy()-fst.get_potential_energy()),
                            abs(trd.get_potential_energy()-fst.get_potential_energy()))
            diff_second=abs(trd.get_potential_energy()-snd.get_potential_energy())
        iters+=1

    #smearing-width convergence test
    diff_primary=100
    diff_second=100
    iters=0
    while diff_primary>rela_tol and diff_second>rela_tol and iters<6: 
        if iters>0:
            sw=sw/2
        atoms=bulk(element)
        calc=GPAW(xc=xc,h=h,kpts=(k,k,k),occupations={'name':'fermi-dirac','width':sw})
        atoms.set_calculator(calc)
        opt.optimize_bulk(atoms,step=0.05,fmax=0.01,location=element+"_"+'bulk'+'/'+'results_sw',extname='{}'.format(sw))
        db_sw.write(atoms,sw=sw)
        if iters>=2:
            fst=db_h.get_atoms(id=iters-1)
            snd=db_h.get_atoms(id=iters)
            trd=db_h.get_atoms(id=iters+1)
            diff_primary=max(abs(snd.get_potential_energy()-fst.get_potential_energy()),
                            abs(trd.get_potential_energy()-fst.get_potential_energy()))
            diff_second=abs(trd.get_potential_energy()-snd.get_potential_energy())
        iters+=1
    
    parprint('converged h = {}'.format(h))
    parprint('converged k = {}'.format(k))
    parprint('converged smearing width = {}'.format(sw))
