from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor
from ase.db import connect
from pymatgen.io.cif import CifWriter
import numpy as np
from collections import Counter 
from itertools import chain 

def sym_all_slab(element,max_ind,layers,vacuum_layer):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgenall=generate_all_slabs(bulk_pym,max_ind,layers,vacuum_layer,
                                lll_reduce=True,center_slab=True,
                                symmetrize=True,in_unit_planes=True)
    print('Miller Index'+'\t'+'Num of Different Shift')
    slab_M=[]
    slab_layer=[]
    for slab in slabgenall:
        slab_M.append(slab.miller_index)
        slab_ase=AseAtomsAdaptor.get_atoms(slab)
        slab_layer.append(len(np.unique(slab_ase.positions[:,2])))
    slab_M_unique = Counter(chain(*slab_M))
    print(slab_M_unique)
    print(list(slab_M_unique.keys()))
    # for key in list(slab_M_unique.keys()):
    #     print(key)
    #     print(slab_M_unique[key])
    #     print(str(key)+'\t'+str(slab_M_unique[key]))



