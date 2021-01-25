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
    print('Miller Index'+'\t'+'Num of Different Shift(s)')
    slab_M=[]
    for slab in slabgenall:
        slab_M.append([slab.miller_index])
    slab_M_unique = Counter(chain(*slab_M))
    for key in list(slab_M_unique.keys()):
        print(str(key)+'\t'+str(slab_M_unique[key]))