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

def surf_creator(element,ind,layers,vacuum_layer,option='slabgen',max_ind=1,unit=True):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    if option=='slabgen':
        slabgen = SlabGenerator(bulk_pym, ind, layers, vacuum_layer,
                            center_slab=True,lll_reduce=True,in_unit_planes=unit)
        slabs=slabgen.get_slabs()
        slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
    elif option=='gen_all':
        slabgenall=generate_all_slabs(bulk_pym,max_ind,layers,vacuum_layer,
                            lll_reduce=True,center_slab=True,
                            symmetrize=True,in_unit_planes=unit)
        slabs=[slab for slab in slabgenall if slab.miller_index==ind]
        slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
    elif option=='ase':
        slabs=surface(bulk_ase,ind,layers=layers,vacuum=vacuum_layer)
        slabs_symmetric=AseAtomsAdaptor.get_structure(slabs)
    if len(slabs_symmetric) == 0:
        print('No symmetric slab found!')
    else:
        print('No.'+'\t'+'Layers'+'\t'+'Angles'+'\t\t\t\tCell Length')
        fig=plt.figure(figsize=(20,10))
        for n,slab in enumerate(slabs_symmetric):
            slab_ase=AseAtomsAdaptor.get_atoms(slab)
            angles=np.round(slab_ase.get_cell_lengths_and_angles()[3:],decimals=4)
            cell_length=np.round(slab_ase.get_cell_lengths_and_angles()[:3],decimals=4)
            print(str(n)+'\t'+str(len(np.unique(slab_ase.positions[:,2])))+'\t'+str(angles)+'\t'+str(cell_length))
            ax=fig.add_subplot(np.ceil(len(slabs_symmetric)/2),2,n+1)
            plot_slab(slab,ax,adsorption_sites=False,decay=0.25,window=1)
            ax.set_title('{}: No. {}'.format(slab.miller_index,n),{'fontsize':10})
            ax.set_xticks([])
            ax.set_yticks([])

def surf_saver(element,ind,layers,vacuum_layer,option):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgen = SlabGenerator(bulk_pym, ind, layers, vacuum_layer,
                            center_slab=True,lll_reduce=True,in_unit_planes=True)
    slabs=slabgen.get_slabs()
    slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
    slab_to_save=slabs_symmetric[option]
    rep_location=element+'/raw_surf'
    if os.path.isdir(rep_location):
        print('WARNING: '+rep_location+' already exists!')
    os.makedirs(rep_location,exist_ok=True)
    surf_location=element+'/raw_surf/'+str(slab_to_save.miller_index)+'.cif'
    if os.path.isdir(surf_location):
        print('WARNING: '+surf_location+' already exists!')
        print('Raw surface saving fail!')
    else:
        CifWriter(slab_to_save).write_file(surf_location)
        print('Raw surface saving complete!')


