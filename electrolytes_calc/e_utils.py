import os
import sys
import pubchempy as pcp

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