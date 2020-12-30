import os
import sys
def create_dir(element,options=['bulk','surf','ads'],
                surf_struc=['100','110','111'],
                optimized_parameters=['h','k','sw'],
                starting_layer=3,
                init_vac=5):

    #create the element dir
    if os.path.isdir(element):
        print("WARNING: {} directory already exists!".format(element))
    else:
        os.makedirs(element,exist_ok=True)
    
    
    #create the bulk dir
    if os.path.isdir(element+'/'+'bulk'):
        print("WARNING: bulk directory already exists!")
    else:
        os.makedirs(element+'/'+'bulk',exist_ok=True)
    for par in optimized_parameters:
        create_bulk_sub_dir(element,par)
    print("bulk directory complete!")

    #create the surf dir

    if os.path.isdir(element+'/'+'surf'):
        print("WARNING: surf directory already exists!")
    else:
        os.makedirs(element+'/'+'surf',exist_ok=True)
    for struc in surf_struc:
        create_surf_sub_dir(element,struc,starting_layer)
        create_surf_vac_dir(element,struc,init_vac)
    print('surf directory complete!')

    #TO-DO create adsorption dir

def create_surf_vac_dir(element,struc,init_vac):
    if os.path.isdir(element+'/'+'surf'+'/'+struc+'/'+'layer_optimized'):
        print('WARNING: '+element+'/'+'surf'+'/'+struc+'/'+'layer_optimized directory already exists!')
    else:
        os.makedirs(element+'/'+'surf'+'/'+struc+'/'+'layer_optimized',exist_ok=True)
    for vac in range(init_vac,init_vac+6):
        if os.path.isdir(element+'/'+'surf'+'/'+struc+'/'+'layer_optimized'+'/'+'vacuum_'+str(vac)):
            print('WARNING: '+element+'/'+'surf'+'/'+struc+'/'+'layer_optimized'+'/'+'vacuum_'+str(vac)+' directory already exists!')
        else:
            os.makedirs(element+'/'+'surf'+'/'+struc+'/'+'layer_optimized'+'/'+'vacuum_'+str(vac),exist_ok=True)

def create_surf_sub_dir(element,struc,starting_layer):
    if os.path.isdir(element+'/'+'surf'+'/'+struc):
        print('WARNING: '+element+'/'+'surf'+'/'+'{} directory already exists!'.format(struc))
    else:
        os.makedirs(element+'/'+'surf'+'/'+struc,exist_ok=True)
    for layer in range(starting_layer,starting_layer+8): #3-10
        if os.path.isdir(element+'/'+'surf'+'/'+struc+'/'+str(layer)+'x1x1'):
            print('WARNING: '+element+'/'+'surf'+'/'+struc+'/'+str(layer)+'x1x1'+' directory already exists!')
        else:
            os.makedirs(element+'/'+'surf'+'/'+struc+'/'+str(layer)+'x1x1',exist_ok=True)
        
def create_bulk_sub_dir(element,par):
    if os.path.isdir(element+'/'+'bulk'+'/'+'results'+'_'+par):
        print('WARNING: '+element+'/'+'bulk'+'/'+'results_{} directory already exists!'.format(par))
    else:
        os.makedirs(element+'/'+'bulk'+'/'+'results'+'_'+par,exist_ok=True)

    if os.path.isdir(element+'/'+'bulk'+'/'+'results'+'_'+par+'/'+'eos_fit'):
        print('WARNING: '+element+'/'+'bulk'+'/'+'results'+'_'+par+'/'+'eos_fit'+'directory already exists!')
    else:
        os.makedirs(element+'/'+'bulk'+'/'+'results'+'_'+par+'/'+'eos_fit',exist_ok=True)