import os
import sys
import pubchempy as pcp
from ase.data.pubchem import available_conformer_search,pubchem_atoms_conformer_search

def pause():
    input('Press <ENTER> to continue...')

def create_big_dir():
    current_dir=os.getcwd()
    os.chdir(current_dir)
    #create the orig_cif_data and final_database dir
    if os.path.isdir('input_xyz'):
        print("WARNING: input_xyz/ directory already exists!")
    else:
        os.makedirs('input_xyz')
    if os.path.isdir('final_database'):
        print("WARNING: final_database/ directory already exists!")
    else:
        os.makedirs('final_database')
    if os.path.isdir('results'):
        print("WARNING: results/ directory already exists!")
    else:
        os.makedirs('results') 

def create_mol_dir(mol_name):
    current_dir=os.getcwd()
    os.chdir(current_dir)
    c=pcp.get_compounds(mol_name,'name')
    cid=str(c[0].cid)
    num_of_conformer=len(available_conformer_search(cid,'cid'))
    cid_i_ls=[]
    if num_of_conformer==1:
        cid_i=cid+'_'+'0'
        os.makedirs('results/'+cid_i,exist_ok=True)
        if os.path.isdir(cid_i):
            print("WARNING: {}(cid={}) directory already exists!".format(mol_name,cid_i))
        else:
            print('{}(cid={}) directory created!'.format(mol_name,cid_i))
        cid_i_ls.append(cid_i)
    else:
        for i in range(1,num_of_conformer+1):
            cid_i=cid+'_'+str(i)
            os.makedirs(cid_i,exist_ok=True)
            if os.path.isdir(cid_i):
                print("WARNING: {}(cid={}) directory already exists!".format(mol_name,cid_i))
            else:
                print('{}(cid={}) directory created!'.format(mol_name,cid_i))
            cid_i_ls.append(cid_i)
    return cid_i_ls

def create_mol_sub_dir(cid,sub_dir):
    current_dir=os.getcwd()
    os.chdir(current_dir+'/'+'results/'+cid)
    for j in sub_dir:
        if os.path.isdir(j):
            print('WARNING: (cid={}) {} directory already exists!'.format(cid,j))
            continue
        else:
            os.makedirs(j)
            print('(cid={}) {} directory created!'.format(cid,j))
    os.chdir(current_dir)

def mol_pubchem_grabber(cid):
    cid_i=cid[0]
    pure_cid=int(cid_i.split('_')[0])
    current_dir=os.getcwd()
    os.chdir(current_dir+'/'+'input_xyz')
    try:
        mol=pubchem_atoms_conformer_search(cid=pure_cid)
    except:
        print("ERROR: Can't find '{}' in PubChem Database.".format(str(cid)))
        sys.exit()
    c=pcp.get_compounds(pure_cid,'cid')
    synonyms_name=(c[0].synonyms)[0]
    mol_name=synonyms_name.lower().replace(' ','-')
    if len(mol)==1:
        mol[0].write(mol_name+'_'+str(pure_cid)+'_'+'0'+'.xyz',format='xyz')
        print("'"+str(pure_cid)+'0'+"'",'input xyz is saved successfully!')
    else:
        for i,mol_i in enumerate(mol):
            mol_i.write(mol_name+'_'+str(pure_cid)+'_'+str(i+1)+'.xyz',format='xyz')
            print("'"+str(pure_cid)+'_'+str(i+1)+"'",'input xyz is saved successfully!')
    os.chdir(current_dir)