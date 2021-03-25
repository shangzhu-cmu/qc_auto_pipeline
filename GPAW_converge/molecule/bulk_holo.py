from gpaw import GPAW,Mixer,Davidson
from ase.build import bulk
from ase.db import connect
import os
import GPAW_converge.molecule.optimizer as opt
from gpaw.eigensolvers import RMMDIIS
from ase.parallel import parprint
import numpy as np
import sys
from ase.io import read, write
from ase.parallel import paropen, parprint, world
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp

def homo_lumo(element,gpaw_calc,relax_xc,
                    init_magmom=0,
                    solver_fmax=0.01,
                    solver_maxstep=0.04):
    calc_dict=gpaw_calc.__dict__['parameters']
    cid=element.split('_')[-2:]
    cid='_'.join(cid)
    parprint(cid)
    rep_location=cid+'/'+'homo-lumo'+'_results_report.txt'
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    with paropen(rep_location,'a') as f:
        parprint('Parameters:',file=f)
        parprint('\t'+'Materials: '+element,file=f)
        parprint('\t'+'xc: '+calc_dict['xc'].split('-')[0],file=f)
        parprint('\t'+'h: '+str(calc_dict['h']),file=f)
        parprint('\t'+'kpts: '+str(calc_dict['kpts']),file=f)
        parprint('\t'+'sw: '+str(calc_dict['occupations']),file=f)
        parprint('\t'+'spin polarized: '+str(calc_dict['spinpol']),file=f)
        if calc_dict['spinpol']:
            parprint('\t'+'magmom: '+str(init_magmom),file=f)
    f.close()
    #connecting to databse
    print(element)
    db_opt=connect('final_database'+'/'+'bulk_'+relax_xc+'.db')
    db_holo=connect('final_database'+'/'+'homo_lumo.db')
    atoms=db_opt.get_atoms(name=element)
    atoms.set_calculator(gpaw_calc)
    #(atoms,cid,XC,fmax=solver_fmax, maxstep=solver_maxstep, replay_traj=None)
    print(atoms.get_potential_energy())
    (homo,lumo)=calc.get_homo_lumo()
    id=db_final.reserve(name=element)
    if id is None:
        id=db_final.get(name=element).id
        db_final.update(id=id,atoms=atoms,name=element,
                        h=calc_dict['h'],sw=calc_dict['occupations']['width'],xc=calc_dict['xc'],spin=calc_dict['spinpol'],
                        kpts=str(','.join(map(str, calc_dict['kpts']))))
    else:
        db_final.write(atoms,id=id,name=element,
                        h=calc_dict['h'],sw=calc_dict['occupations']['width'],xc=calc_dict['xc'],spin=calc_dict['spinpol'],
                        kpts=str(','.join(map(str, calc_dict['kpts']))))

def mol_builder(element):
    location='input_xyz'+'/'+element+'.xyz'
    atoms=read(location)
    pos = atoms.get_positions()
    xl = max(pos[:,0])-min(pos[:,0])+15
    yl = max(pos[:,1])-min(pos[:,1])+15
    zl = max(pos[:,2])-min(pos[:,2])+15
    maxlength=max([xl,yl,zl])
    atoms.set_cell((maxlength,maxlength,maxlength))
    atoms.center()
    atoms.set_pbc([False,False,False])
    return atoms