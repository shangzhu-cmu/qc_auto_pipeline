from ase.optimize import BFGS
import numpy as np
from ase.dft.bee import BEEFEnsemble
from ase.parallel import parprint

def relax_single(atoms,name,sub_dir,fmax=0.01, maxstep=0.04):
    file_dir_name='results/'+name+'/'+sub_dir+'/'+'mol'
    atoms.calc.set(txt=file_dir_name+'.txt')
    dyn=BFGS(atoms=atoms,trajectory=file_dir_name+'.traj',logfile = file_dir_name+'.log',restart=file_dir_name+'qn.pckl',maxstep=maxstep)
    dyn.run(fmax=fmax)
    #atoms.calc.write(file_dir_name+'.gpw')
    #Writing ensemble energies to file 
    if atoms.calc.get_xc_functional()=='BEEF-vdW':
        ens = BEEFEnsemble(atoms.calc)
        ens_material = ens.get_ensemble_energies()
        with open(file_dir_name+'_Ensemble_Energies.txt','w+') as file:
            file.write(file_dir_name+'\n')
            for i in range(len(ens_material)):
                file.write(str(ens_material[i])+'\n')
    return file_dir_name
    

def SPE_calc(atoms,name,save_gpw=True):
    file_path='results/'+name
    atoms.calc.set(txt=file_path+'.txt')
    atoms.get_potential_energy()
    if save_gpw == True:
        atoms.calc.write(file_path+'.gpw',mode='all')
    return file_path
