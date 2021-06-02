from gpaw import GPAW
from ase.build import molecule,bulk
from ase import Atoms
import GPAW_converge.molecule.optimizer as opt
from ase.db import connect

def ele_calc(atoms,element,gpaw_calc,init_mag=0,uc_size=15,solver_fmax=0.01,solver_maxstep=0.04):
    calc_dict=gpaw_calc.__dict__['parameters']
    XC=calc_dict['xc'].split('-')[0]
    atoms.set_calculator(gpaw_calc)
    atoms.set_initial_magnetic_moments(init_mag)
    opt.relax_single(atoms,element,XC,fmax=solver_fmax, maxstep=solver_maxstep, replay_traj=None)
    db_element=connect('final_database'+'/'+'element_'+calc_dict['xc']+'.db')
    id=db_element.reserve(name=element)
    if id is None:
        id=db_element.get(name=element).id
        db_element.update(id=id,atoms=atoms,name=element,
                        h=calc_dict['h'],sw=calc_dict['occupations']['width'],xc=calc_dict['xc'],spin=calc_dict['spinpol'],
                        kpts=str(','.join(map(str, calc_dict['kpts']))))
    else:
        db_element.write(atoms,id=id,name=element,
                        h=calc_dict['h'],sw=calc_dict['occupations']['width'],xc=calc_dict['xc'],spin=calc_dict['spinpol'],
                        kpts=str(','.join(map(str, calc_dict['kpts']))))


def mol_builder(element,uc_size):
    mol=molecule(element)
    pos = mol.get_positions()
    xl = max(pos[:,0])-min(pos[:,0])+uc_size
    yl = max(pos[:,1])-min(pos[:,1])+uc_size
    zl = max(pos[:,2])-min(pos[:,2])+uc_size
    maxlength=max([xl,yl,zl])
    mol.set_cell((maxlength,maxlength,maxlength))
    mol.center()
    mol.set_pbc([False,False,False])
    return mol

def bulk_builder(element,uc_size):
    atoms=Atoms(element)
    pos = atoms.get_positions()
    xl = max(pos[:,0])-min(pos[:,0])+uc_size
    yl = max(pos[:,1])-min(pos[:,1])+uc_size
    zl = max(pos[:,2])-min(pos[:,2])+uc_size
    maxlength=max([xl,yl,zl])
    atoms.set_cell((maxlength,maxlength,maxlength))
    atoms.center()
    atoms.set_pbc([False,False,False])
    return atoms

def atoms_builder(element,uc_size):
    atoms=Atoms(element)
    pos = atoms.get_positions()
    xl = max(pos[:,0])-min(pos[:,0])+uc_size
    yl = max(pos[:,1])-min(pos[:,1])+uc_size
    zl = max(pos[:,2])-min(pos[:,2])+uc_size
    maxlength=max([xl,yl,zl])
    atoms.set_cell((maxlength,maxlength,maxlength))
    atoms.center()
    atoms.set_pbc([False,False,False])
    return atoms


