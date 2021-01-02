from gpaw import GPAW
from ase.build import bulk
from ase.io import write
from ase.db import connect
import os
import optimizer as opt
from ase.parallel import parprint
import numpy as np
import sys
from ase.io import read

def bulk_auto_conv(element,a0=None,struc=None,h=0.16,k=6,xc='PBE',sw=0.1,rela_tol=10*10**(-3),cif=False,temp_print=True):
    db_h=connect(element+"/"+'bulk'+'/'+'grid_converge.db')
    db_k=connect(element+"/"+'bulk'+'/'+'kpts_converge.db')
    db_sw=connect(element+"/"+'bulk'+'/'+'sw_converge.db')
    diff_primary=100
    diff_second=100
    grid_iters=0
    #start with grid spacing convergence
    h_ls=[]
    while (diff_primary>rela_tol or diff_second>rela_tol) and grid_iters <= 6:
        if grid_iters>0:
            h=np.round(h-0.02,decimals=2)
        atoms=bulk_builder(element,cif,struc,a0)
        calc=GPAW(xc=xc,h=h,kpts=(k,k,k),occupations={'name':'fermi-dirac','width':sw})
        atoms.set_calculator(calc)
        opt.optimize_bulk(atoms,step=0.05,fmax=0.01,location=element+"/"+'bulk'+'/'+'results_h',extname='{}'.format(h))
        db_h.write(atoms,h=h)
        if grid_iters>=2:
            fst=db_h.get_atoms(id=grid_iters-1)
            snd=db_h.get_atoms(id=grid_iters)
            trd=db_h.get_atoms(id=grid_iters+1)
            diff_primary=max(abs(snd.get_potential_energy()-fst.get_potential_energy()),
                            abs(trd.get_potential_energy()-fst.get_potential_energy()))
            diff_second=abs(trd.get_potential_energy()-snd.get_potential_energy())
            temp_output_printer(db_h,grid_iters,'h',temp_print)
        h_ls.append(h)
        grid_iters+=1
    if grid_iters>=6:
        if diff_primary>rela_tol or diff_second>rela_tol:
            parprint("WARNING: Max GRID iterations reached! System may not converged.")
            parprint("Possible Error: Incorrect Lattice Parameters, Inappropriate Starting Grid Size.")
            parprint("Computation Suspended!")
            sys.exit()
    h=h_ls[-3]
    #kpts convergence
    diff_primary=100
    diff_second=100
    k_iters=0
    k_ls=[]
    while (diff_primary>rela_tol or diff_second>rela_tol) and k_iters <= 6: 
        if k_iters>0:
            k=int(k+2)
        atoms=bulk_builder(element,cif,struc,a0)
        #atoms=bulk(element,struc,a=a0)
        calc=GPAW(xc=xc,h=h,kpts=(k,k,k),occupations={'name':'fermi-dirac','width':sw})
        atoms.set_calculator(calc)
        opt.optimize_bulk(atoms,step=0.05,fmax=0.01,location=element+"/"+'bulk'+'/'+'results_k',extname='{}'.format(k))
        db_k.write(atoms,kpts=k)
        if k_iters>=2:
            fst=db_k.get_atoms(id=k_iters-1)
            snd=db_k.get_atoms(id=k_iters)
            trd=db_k.get_atoms(id=k_iters+1)
            diff_primary=max(abs(snd.get_potential_energy()-fst.get_potential_energy()),
                            abs(trd.get_potential_energy()-fst.get_potential_energy()))
            diff_second=abs(trd.get_potential_energy()-snd.get_potential_energy())
            temp_output_printer(db_k,k_iters,'kpts',temp_print)
        k_iters+=1
        k_ls.append(k)
    if k_iters>=6:
        if diff_primary>rela_tol or diff_second>rela_tol:
            parprint("WARNING: Max Kpts iterations reached! System may not converged.")
            #parprint("Possible Error: Incorrect Lattice Parameters, Inappropriate Starting Kpts Size.")
            parprint("Computation Suspended!")
            sys.exit()
    k=k_ls[-3]
    #smearing-width convergence test
    diff_primary=100
    diff_second=100
    sw_iters=0
    sw_ls=[]
    while (diff_primary>rela_tol or diff_second>rela_tol) and sw_iters <= 6: 
        if sw_iters>0:
            sw=sw/2
        atoms=bulk_builder(element,cif,struc,a0)
        #atoms=bulk(element,struc,a=a0)
        calc=GPAW(xc=xc,h=h,kpts=(k,k,k),occupations={'name':'fermi-dirac','width':sw})
        atoms.set_calculator(calc)
        opt.optimize_bulk(atoms,step=0.05,fmax=0.01,location=element+"/"+'bulk'+'/'+'results_sw',extname='{}'.format(sw))
        db_sw.write(atoms,sw=sw)
        if sw_iters>=2:
            fst=db_sw.get_atoms(id=sw_iters-1)
            snd=db_sw.get_atoms(id=sw_iters)
            trd=db_sw.get_atoms(id=sw_iters+1)
            diff_primary=max(abs(snd.get_potential_energy()-fst.get_potential_energy()),
                            abs(trd.get_potential_energy()-fst.get_potential_energy()))
            diff_second=abs(trd.get_potential_energy()-snd.get_potential_energy())
            temp_output_printer(db_sw,sw_iters,'sw',temp_print)
        sw_iters+=1
        sw_ls.append(sw)
    if sw_iters>=6:
        if diff_primary>rela_tol or diff_second>rela_tol:
            parprint("WARNING: Max Smearing-Width iterations reached! System may not converged.")
            #parprint("Possible Error: Incorrect Lattice Parameters, Inappropriate Starting Kpts Size.")
            parprint("Computation Suspended!")
            sys.exit()
    sw=sw_ls[-3]
    final_atom=db_sw.get_atoms(id=len(db_sw)-2)
    parprint('converged h = {}'.format(h))
    parprint('converged k = {}'.format(k))
    parprint('converged smearing width = {}'.format(sw))
    parprint('converged final lattice constant a = {} Ang'.format(np.round(final_atom.cell[0][1]*2,decimals=5))) ## TO-DO: NEED a comprehensive output
    parprint('converged final potential energy e = {} eV'.format(np.round(final_atom.get_potential_energy(),decimals=5)))

def bulk_builder(element,cif,struc,a0):
    if cif == False:
        atoms=bulk(element,struc,a=a0)
    else:
        location='orig_cif_data'+'/'+element+'.cif'
        atoms=read(location)
    return atoms

def temp_output_printer(db,iters,key,option=False):
    fst_r=db.get(iters-1)
    snd_r=db.get(iters)
    trd_r=db.get(iters+1)
    if option==True:
        parprint('2nd{}-1st{}'.format(snd_r[key],fst_r[key]),
                '=',
                np.round(abs(snd_r['energy']-fst_r['energy']),decimals=5),'eV')
        parprint('3rd{}-1st{}'.format(trd_r[key],fst_r[key]),
                '=',
                np.round(abs(trd_r['energy']-fst_r['energy']),decimals=5),'eV')
        parprint('3rd{}-2nd{}'.format(trd_r[key],snd_r[key]),
                '=',
                np.round(abs(trd_r['energy']-snd_r['energy']),decimals=5),'eV')