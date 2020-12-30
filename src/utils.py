import os
import sys
def create_dir(element,options=['bulk','surf','ads'],
                surf_struc=['100','110','111'],
                optimized_parameters=['h','k','sw'],
                starting_layer=3,
                init_vac=5):

    #create the element dir
    try:
        os.makedirs(element)
    except:
        print("ERROR: {} directory already exists!".format(element))
        print("Creation Suspended")
        sys.exit()
    
    #create the bulk dir
    try:
        os.makedirs(element+"/"+'bulk')
    except:
        print("ERROR: bulk directory already exists!")
        print("Creation Suspended")
        sys.exit()
    for par in optimized_parameters:
        create_bulk_sub_dir(element,par)
    print("bulk directory complete!")

    #create the surf dir
    try:
        os.makedirs(element+"/"+'surf')
    except:
        print("ERROR: surf directory already exists!")
        print("Creation Suspended")
        sys.exit()
    for struc in surf_struc:
        create_surf_sub_dir(element,struc,starting_layer))
        create_surf_vac_dir(element,struc,init_vac)
    print("surf directory complete!")

    #TO-DO create adsorption dir

def create_surf_vac_dir(element,struc,init_vac):
    try:
        os.makedirs(element+"/"+'surf'+'/'+struc+'/'+'layer_optimized')
    except:
        print("ERROR: "+element+"/"+'surf'+'/'+struc+'/'+'layer_optimized directory already exists!')
        print("Creation Suspended")
        sys.exit()
    for str(vac) in range(init_vac,init_vac+6):
        try:
            os.makedirs(element+"/"+'surf'+'/'+struc+'/'+'layer_optimized'+'/'+'vacuum_'+vac)
        except:
            print("ERROR: "+element+"/"+'surf'+'/'+struc+'/'+'layer_optimized'+'/'+'vacuum_'+vac+" directory already exists!")
            print("Creation Suspended")
            sys.exit()   

def create_surf_sub_dir(element,struc,starting_layer):
    try:
        os.makedirs(element+"/"+'surf'+'/'+struc)
    except:
        print("ERROR: "+element+"/"+"surf"+"/"+"{} directory already exists!".format(struc))
        print("Creation Suspended")
        sys.exit()
    for str(layer) in range(starting_layer,starting_layer+8): #3-10
        try:
            os.makedirs(element+"/"+'surf'+'/'+struc+'/'+layer+"x1x1")
        except:
            print("ERROR: "+element+"/"+'surf'+'/'+struc+'/'+layer+"x1x1"+" directory already exists!")
            print("Creation Suspended")
            sys.exit()
        
def create_bulk_sub_dir(element,par):
    try:
        os.makedirs(element+"/"+'bulk'+'/'+'results'+"_"+par)
    except:
        print("ERROR: "+element+"/"+"bulk"+"/"+"results_{} directory already exists!".format(par))
        print("Creation Suspended")
        sys.exit()
    try:
        os.makedirs(element+"/"+'bulk'+'/'+'results'+"_"+par+'/'+'eos_fit')
    except:
        print("ERROR: "+element+"/"+'bulk'+'/'+'results'+"_"+par+'/'+'eos_fit'+"directory already exists!")
        print("Creation Suspended")
        sys.exit()