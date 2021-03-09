import os
import sys
import pubchempy as pcp
from ase.data.pubchem import pubchem_atoms_search

def pause():
    input('Press <ENTER> to continue...')

def create_mol_dir(mol_name,sub_dir=['PBE','BEEF']):
    current_dir=os.getcwd()
    os.chdir(current_dir)
    c=pcp.get_compounds(mol_name,'name')
    cid=str(c[0].cid)
    #create the orig_cif_data and final_database dir
    if os.path.isdir(cid):
        print("WARNING: {}(cid={}) directory already exists!".format(mol_name,cid))
        sys.exit()
    else:
        os.makedirs(cid)
        for i in sub_dir:
            os.makedirs(cid+'/'+i)
        return int(cid) 


def create_big_dir():
    current_dir=os.getcwd()
    os.chdir(current_dir)
    #create the orig_cif_data and final_database dir
    if os.path.isdir('input_xyz'):
        print("WARNING: input_xyz directory already exists!")
        sys.exit()
    else:
        os.makedirs('input_xyz')
    if os.path.isdir('final_database'):
        print("WARNING: final_database directory already exists!")
        sys.exit()
    else:
        os.makedirs('final_database')

def mol_pubchem_grabber(cid):
    try:
        mol=pubchem_atoms_search(cid=cid,silent=True)
        c=pcp.get_compounds(cid,'cid')
        synonyms_name=(c[0].synonyms)[0]
        mol_name=synonyms_name.lower().replace(' ','-')
        mol[0].write('./input_xyz/'+mol_name+'_'+str(cid)+'.xyz')
        print("'"+mol_name+"'",'input xyz is saved successfully!')
    except:
        print("ERROR: Can't find '{}' in PubChem Database.".format(str(cid)))
