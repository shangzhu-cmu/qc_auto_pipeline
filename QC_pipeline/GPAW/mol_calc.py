from gpaw import restart
from ase.db import connect
import QC_pipeline.GPAW.optimizer as opt
import numpy as np
from ase.parallel import paropen
from ase.io import read
class GPAW_mol_calculator:
    def __init__(self,element):
        self.element=element
        
    def relax_mol(self,
                calculator,
                sub_dir=None,
                init_magmom=0,
                solver_fmax=0.01,
                solver_maxstep=0.04,
                ):

        calc_dict=calculator.__dict__['parameters']
        cid=self.element.split('_')[-2:]
        cid='_'.join(cid)

        if sub_dir is None:
            sub_dir=calc_dict['xc']
        
        self.atoms.set_calculator(calculator)
        if calc_dict['spinpol'] == True:
            self.atoms.set_initial_magnetic_moments(init_magmom*np.ones(len(self.atoms)))

        opt.relax_single(self.atoms,cid,sub_dir,solver_fmax,solver_maxstep)
        self.file_dir_name='results/'+self.element+'/'+sub_dir+'/'+'mol'
        self.database_save('relaxed_'+sub_dir)
        # return self.atoms

    def homo_lumo_calc(self,
                    relax_sub_dir,
                    calculator=None,
                    file_name='mol',
                    mode='occupied',#TWO OTHER MODE: "add_bands", "unoccupied"
                    add_bands=15, 
                    above_lumo_percent=0.3, #percentage of the range above lumo
                    ):
        cid=self.element.split('_')[-2:]
        cid='_'.join(cid)
        if mode == 'occupied':
            file_prev='results/'+cid+'/'+relax_sub_dir+'/'+'mol'
            self.atoms = restart(file_prev+'.gpw')[0]
            if calculator == None:
                raise RuntimeError('No HOMO LUMO Calculator.')
            self.atoms.set_calculator(calculator)
            opt.SPE_calc(self.atoms,name=cid+'/'+'homo-lumo'+'/'+file_name+'_occupied')
        elif mode == 'add_bands':
            file_prev='results/'+cid+'/'+'homo-lumo'+'/'+file_name+'_occupied'
            nbands=nbands_finder(file_prev+'.txt')
            self.atoms, calculator = restart(file_prev+'.gpw',nbands=nbands+add_bands)
            self.atoms.set_calculator(calculator)
            opt.SPE_calc(self.atoms,name=cid+'/'+'homo-lumo'+'/'+file_name+'_add_bands')
        elif mode == 'unoccupied':
            file_prev='results/'+cid+'/'+'homo-lumo'+'/'+file_name+'_add_bands'
            eigen_arr=aboveLUMO_finder(file_prev+'.txt')
            aboveLUMO=np.abs(max(eigen_arr)-min(eigen_arr))*above_lumo_percent
            self.atoms, calculator = restart(file_prev+'.gpw')
            calculator['parameters']['convergence']['bands']='CBM+'+str(aboveLUMO)
            self.atoms.set_calculator(calculator)
            opt.SPE_calc(self.atoms,name=cid+'/'+'homo-lumo'+'/'+file_name+'_unoccupied')
            self.database_save('HOLO_'+file_name)
        else:
            raise NameError('mode not definied. Available modes: occupied, add_bands, unoccupied')
        
    def database_save(self,name):
        
        db_final=connect('final_database'+'/'+name+'.db')
        id=db_final.reserve(name=self.element)
        if id is None:
            id=db_final.get(name=self.element).id
            db_final.update(id=id,atoms=self.atoms,name=self.element,
                            gpw_dir=self.file_dir_name+'.gpw')
        else:
            db_final.write(self.atoms,id=id,name=self.element,
                            gpw_dir=self.file_dir_name+'.gpw')

    def bulk_builder(self,box_size):
        location='input_xyz'+'/'+self.element+'.xyz'
        self.atoms=read(location)
        pos = self.atoms.get_positions()
        xl = max(pos[:,0])-min(pos[:,0])+box_size
        yl = max(pos[:,1])-min(pos[:,1])+box_size
        zl = max(pos[:,2])-min(pos[:,2])+box_size
        maxlength=max([xl,yl,zl])
        self.atoms.set_cell((maxlength,maxlength,maxlength))
        self.atoms.center()
        self.atoms.set_pbc([False,False,False])
        return self.atoms

def nbands_finder(file_name):
    file = paropen(file_name,'r')
    for line in file:
        line = line.split()
        if not line:
            continue
        if len(line) >= 5:
            if line[0] == 'Number' and line[1] == 'of' and line[2] == 'bands' and line[3] == 'in' and line[4] == 'calculation:':
                return int(line[-1])

def aboveLUMO_finder(file_name):
    file = paropen(file_name,'r')
    eigen=[]
    for line in file:
        line = line.split()
        if not line:
            continue
        if len(line) >= 3:
            if line[2] == '0.00000':
                eigen.append(float(line[1]))
    return np.array(eigen)