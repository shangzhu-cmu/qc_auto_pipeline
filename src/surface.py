from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor
from ase.db import connect
from ase.build import surface
from pymatgen.io.cif import CifWriter
import numpy as np
from collections import Counter 
from itertools import chain 
from pymatgen.analysis.adsorption import plot_slab
from matplotlib import pyplot as plt
from ase.visualize.plot import plot_atoms
import os

def sym_all_slab(element,max_ind,layers,vacuum_layer):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgenall=generate_all_slabs(bulk_pym,max_ind,layers,vacuum_layer,
                                lll_reduce=True,center_slab=True,
                                symmetrize=True,in_unit_planes=True)
    print('Miller Index'+'\t'+'Num of Different Shift(s)')
    slab_M=[]
    for slab in slabgenall:
        slab_M.append([slab.miller_index])
    slab_M_unique = Counter(chain(*slab_M))
    for key in list(slab_M_unique.keys()):
        print(str(key)+'\t'+str(slab_M_unique[key]))

def surf_creator(element,ind,layers,vacuum_layer,option='slabgen',max_ind=1,unit=True,order=0,save=False):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    if option=='slabgen':
        slabgen = SlabGenerator(bulk_pym, ind, layers, vacuum_layer,
                            center_slab=True,lll_reduce=True,in_unit_planes=unit)
        slabs=slabgen.get_slabs()
        slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
        if len(slabs_symmetric) == 0:
            print('No symmetric slab found!')
        else:
            print('No.'+'\t'+'Layers'+'\t'+'Angles'+'\t\t\t\tCell Length')
            fig=plt.figure(figsize=(8,8))
            for n,slab in enumerate(slabs_symmetric):
                slab_ase=AseAtomsAdaptor.get_atoms(slab)
                angles=np.round(slab_ase.get_cell_lengths_and_angles()[3:],decimals=4)
                cell_length=np.round(slab_ase.get_cell_lengths_and_angles()[:3],decimals=4)
                print(str(n)+'\t'+str(len(np.unique(slab_ase.positions[:,2])))+'\t'+str(angles)+'\t'+str(cell_length))
                ax=fig.add_subplot(np.ceil(len(slabs_symmetric)/2),2,n+1)
                plot_slab(slab,ax,adsorption_sites=False,decay=0.25,window=1)
                ax.set_title('{}: No. {}'.format(slab.miller_index,n),{'fontsize':20})
                ax.set_xticks([])
                ax.set_yticks([])
        if save:
            surf_saver(element,slabs_symmetric[order],ind)
    elif option=='ase':
        slab_ase=surface(bulk_ase,ind,layers=layers,vacuum=vacuum_layer)
        print('No.'+'\t'+'Layers'+'\t'+'Angles'+'\t\t\t\tCell Length')
        angles=np.round(slab_ase.get_cell_lengths_and_angles()[3:],decimals=4)
        cell_length=np.round(slab_ase.get_cell_lengths_and_angles()[:3],decimals=4)
        print(str(0)+'\t'+str(len(np.unique(slab_ase.positions[:,2])))+'\t'+str(angles)+'\t'+str(cell_length))
        fig=plt.figure(figsize=(8,8))
        ax=fig.add_subplot(111)
        plot_atoms(slab_ase,ax=ax)
        ax.set_title('ASE created: {}'.format(str(ind)),{'fontsize':20})
        ax.set_xticks([])
        ax.set_yticks([])
        if save:
            slab_struc=AseAtomsAdaptor.get_structure(slab_ase)
            surf_saver(element,slab_struc,ind)

def surf_saver(element,slab_to_save,ind):
    rep_location=element+'/raw_surf'
    if os.path.isdir(rep_location):
        print('WARNING: '+rep_location+' already exists!')
    os.makedirs(rep_location,exist_ok=True)
    surf_location=element+'/raw_surf/'+str(ind)+'.cif'
    if os.path.isdir(surf_location):
        print('WARNING: '+surf_location+' already exists!')
        print('Raw surface saving fail!')
    else:
        CifWriter(slab_to_save).write_file(surf_location)
        print('Raw surface saving complete!')


