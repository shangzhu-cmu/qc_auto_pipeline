from gpaw import GPAW,Mixer,Davidson
from ase.build import bulk
from ase.db import connect
import os
import GPAW_converge.molecule.optimizer as opt
from ase.parallel import parprint
import numpy as np
import sys
from ase.io import read, write
from ase.parallel import paropen, parprint, world
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp

def bulk_energy(element,gpaw_calc,
                    init_magmom=0,
                    solver_fmax=0.01,
                    solver_maxstep=0.04):
    calc_dict=gpaw_calc.__dict__['parameters']
    cid=element.split('_')[-2:]
    orig_atom=bulk_builder(element)
    XC=calc_dict['xc'].split('-')[0]
    rep_location=(cid+'/'+XC+'/'+'results_report.txt')
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    with paropen(rep_location,'a') as f:
        parprint('Parameters:',file=f)
        parprint('\t'+'Materials: '+element,file=f)
        parprint('\t'+'Lattice constants: '+str(np.round(orig_atom.get_cell_lengths_and_angles()[:3],decimals=5))+'Ang',file=f)    
        parprint('\t'+'Lattice angles: '+str(np.round(orig_atom.get_cell_lengths_and_angles()[3:],decimals=5))+'Degree',file=f)
        parprint('\t'+'xc: '+XC,file=f)
        parprint('\t'+'h: '+str(calc_dict['h']),file=f)
        parprint('\t'+'kpts: '+str(calc_dict['kpts']),file=f)
        parprint('\t'+'sw: '+str(calc_dict['occupations']),file=f)
        parprint('\t'+'spin polarized: '+str(calc_dict['spinpol']),file=f)
        if calc_dict['spinpol']:
            parprint('\t'+'magmom: '+str(init_magmom),file=f)
    f.close()
    #connecting to databse
    db_final=connect('final_database'+'/'+'bulk_'+calc_dict['xc']+'.db')
    atoms=bulk_builder(element)
    atoms.set_calculator(gpaw_calc)
    opt.relax(atoms,cid,XC,fmax=solver_fmax, maxstep=solver_maxstep, replay_traj=None)
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


def bulk_builder(element):
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

def temp_output_printer(db,iters,key,location):
    fst_r=db.get(iters-1)
    snd_r=db.get(iters)
    trd_r=db.get(iters+1)
    with paropen(location,'a') as f:
        parprint('Optimizing parameter: '+key,file=f)
        parprint('\t'+'1st: '+str(fst_r[key])+' 2nd: '+str(snd_r[key])+' 3rd: '+str(trd_r[key]),file=f)
        parprint('\t'+'2nd-1st: '+str(np.round(abs(snd_r['energy']-fst_r['energy']),decimals=5))+'eV',file=f)
        parprint('\t'+'3rd-1st: '+str(np.round(abs(trd_r['energy']-fst_r['energy']),decimals=5))+'eV',file=f)
        parprint('\t'+'3rd-2nd: '+str(np.round(abs(trd_r['energy']-snd_r['energy']),decimals=5))+'eV',file=f)
    f.close()

def mp2kdens(atoms,kpts):
    recipcell=atoms.get_reciprocal_cell()
    kptdensity_ls=[]
    for i in range(len(kpts)):
        kptdensity = kpts[i]/(2 * np.pi * np.sqrt((recipcell[i]**2).sum()))
        kptdensity_ls.append(np.round(kptdensity,decimals=4))
    return kptdensity_ls