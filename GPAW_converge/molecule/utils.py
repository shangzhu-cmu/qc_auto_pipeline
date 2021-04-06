import os
import sys
import pubchempy as pcp
from ase.data.pubchem import pubchem_atoms_search,available_conformer_search,pubchem_atoms_conformer_search
from ase.io import write
def pause():
    input('Press <ENTER> to continue...')

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

def create_mol_dir(mol_name):
    current_dir=os.getcwd()
    os.chdir(current_dir)
    c=pcp.get_compounds(mol_name,'name')
    cid=str(c[0].cid)
    num_of_conformer=len(available_conformer_search(cid,'cid'))
    cid_i_ls=[]
    if num_of_conformer==1:
        cid_i=cid+'_'+'0'
        if os.path.isdir(cid_i):
            print("WARNING: {}(cid={}) directory already exists!".format(mol_name,cid_i))
            sys.exit()
        else:
            os.makedirs(cid_i)
            print('{}(cid={}) directory created!'.format(mol_name,cid_i))
            cid_i_ls.append(cid_i)
    else:
        for i in range(1,num_of_conformer+1):
            cid_i=cid+'_'+str(i)
            if os.path.isdir(cid_i):
                print("WARNING: {}(cid={}) directory already exists!".format(mol_name,cid_i))
                continue
            else:
                os.makedirs(cid_i)
                print('{}(cid={}) directory created!'.format(mol_name,cid_i))
                cid_i_ls.append(cid_i)
    return cid_i_ls

def create_xc_dir(cid,sub_dir=['PBE,BEEF']):
    current_dir=os.getcwd()
    os.chdir(current_dir+'/'+cid)
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
    for mol_i in mol:
        mol_i.write(mol_name+'_'+str(cid)+'.xyz',format='xyz')
        print("'"+str(cid)+"'",'input xyz is saved successfully!')
    os.chdir(current_dir)


def create_converg_dir(mol_name,sub_dir=['PBE','BEEF'],convergence=False,parameters=['h','k','sw']):
    current_dir=os.getcwd()
    os.chdir(current_dir)
    c=pcp.get_compounds(mol_name,'name')
    cid_ls=[str(i) for i in c.cid]
    for i in cid_ls:
    #create the orig_cif_data and final_database dir
        if os.path.isdir(cid):
            print("WARNING: {}(cid={}) directory already exists!".format(mol_name,cid))
            sys.exit()
        else:
            os.makedirs(cid)
            for i in sub_dir:
                os.makedirs(cid+'/'+i)
                if convergence:
                    for j in parameters:
                        os.makedirs(cid+'/'+i+'/results_'+j)
            return int(cid) 