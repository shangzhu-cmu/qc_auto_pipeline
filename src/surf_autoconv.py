from gpaw import GPAW
from ase.constraints import FixAtoms
from ase.build import surface
from ase.io import write,read
from ase.db import connect
import os
import optimizer as opt
from ase.parallel import parprint,paropen,world
import numpy as np
import re
import sys
import time
import copy as cp
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
from shutil import copyfile
import glob
from gpaw import Davidson
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from gpaw import Mixer


def surf_auto_conv(element,struc,init_layer=5,vac=5,fix_layer=2,rela_tol=5,temp_print=True,generator='pymatgen'):
    #convert str ind to tuple
    m_ind=tuple(map(int,struc))

    #create report
    rep_location=(element+'/'+'surf'+'/'+struc+'_results_report.txt')
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    
    #check the optimized bulk object
    if not os.path.isfile('final_database/bulk.db'):
        with paropen(rep_location,'a') as f:
            parprint('ERROR: bulk database has not been established!',file=f)
            parprint('Surface Convergence Computation Suspended!',file=f)
        f.close()
        sys.exit()
    else:
        db_bulk=connect('final_database'+'/'+'bulk.db')
        try:
            opt_bulk=db_bulk.get_atoms(name=element)
        except:
            with paropen(rep_location,'a') as f:
                parprint('ERROR: No Optimized Bulk Object Found!',file=f)
                parprint('Surface Convergence Computation Suspended!',file=f)
            f.close()
            sys.exit()
    
    #get the optimized bulk object and converged parameters
    pymatgen_bulk=AseAtomsAdaptor.get_structure(opt_bulk)
    xc=db_bulk.get(name=element).xc
    h=db_bulk.get(name=element).h
    k_density=db_bulk.get(name=element).k_density
    kpts=[int(i) for i in (db_bulk.get(name=element).kpts).split(',')]
    sw=db_bulk.get(name=element).sw
    
    #print out parameters
    with paropen(rep_location,'a') as f:
        parprint('Initial Parameters:',file=f)
        parprint('\t'+'Materials: '+element,file=f)
        parprint('\t'+'Miller Index: '+str(m_ind),file=f)
        parprint('\t'+'Actual Layer: '+str(init_layer),file=f)
        parprint('\t'+'Vacuum length: '+str(vac)+'Ang',file=f)
        parprint('\t'+'Fixed layer: '+str(fix_layer),file=f)
        parprint('\t'+'xc: '+xc,file=f)
        parprint('\t'+'h: '+str(h),file=f)
        parprint('\t'+'k_density: '+str(k_density),file=f)
        parprint('\t'+'kpts: '+str(kpts),file=f)
        parprint('\t'+'sw: '+str(sw),file=f)
        parprint('\t'+'rela_tol: '+str(rela_tol)+'%',file=f)
    f.close()

    #optimize the layers
    ##connect to the layer convergence database
    db_layer=connect(element+'/'+'surf'+'/'+struc+'/'+'layer_converge.db')
    diff_primary=100
    diff_second=100
    iters=len(db_layer)-1
    act_layer_ls=[]
    sim_layer_ls=[]
    sim_layer=1
    if generator=='pymatgen':
        slabgen = SlabGenerator(pymatgen_bulk, m_ind, sim_layer, sim_layer*2, center_slab=True, lll_reduce=True, in_unit_planes=True)
        slabs=slabgen.get_slabs() #this only take the first structure
        slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
        slab=AseAtomsAdaptor.get_atoms(slabs_symmetric[0]) #convert to ase structure
    elif generator=='ase':
        slab=surface(opt_bulk,m_ind,layers=sim_layer,vacuum=vac)
    elif generator=='import':
        slab=read(element+'/raw_surf/'+str(m_ind)+'_'+str(init_layer)+'.cif')
    actual_layer=len(np.unique(np.round(slab.positions[:,2],decimals=4)))
    while (diff_primary>rela_tol or diff_second>rela_tol) and iters <= 5:
        while actual_layer != init_layer:
            sim_layer+=1
            if generator=='pymatgen':
                slabgen = SlabGenerator(pymatgen_bulk, m_ind, sim_layer, sim_layer*2, center_slab=True, lll_reduce=True, in_unit_planes=True)
                slabs=slabgen.get_slabs() #this only take the first structure
                slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
                slab=AseAtomsAdaptor.get_atoms(slabs_symmetric[0]) #convert to ase structure
            elif generator=='ase':
                slab=surface(opt_bulk,m_ind,layers=sim_layer,vacuum=vac)
            else:
                with paropen(rep_location,'a') as f:
                    parprint('ERROR: The number of layers of the imported surface is not correct!',file=f)
                    parprint('\t'+'Actual Layer: '+str(actual_layer),file=f)
                    parprint('\t'+'Desired Layer: '+str(init_layer),file=f)
                    parprint('Computation Suspended!',file=f)
                    sys.exit()    
            actual_layer=len(np.unique(np.round(slab.positions[:,2],decimals=4)))
            if actual_layer > init_layer:
                with paropen(rep_location,'a') as f:
                    parprint('ERROR: Actual number of layers is greater than the desired number of layers.',file=f)
                    parprint('\t'+'Actual Layer: '+str(actual_layer),file=f)
                    parprint('\t'+'Desired Layer: '+str(init_layer),file=f)
                    parprint('Computation Suspended!',file=f)
                    sys.exit()    
        current_vac=slab.cell.lengths()[-1]-max(slab.positions[:,2])
        if current_vac != vac:
            slab.center(vacuum=vac,axis=2)
        fix_mask=np.round(slab.positions[:,2],decimals=4) <= np.unique(np.round(slab.positions[:,2],decimals=4))[fix_layer-1]
        slab.set_constraint(FixAtoms(mask=fix_mask))
        slab.set_pbc([1,1,0])
        kpts=kdens2mp(slab,kptdensity=k_density,even=True)
        slab_length=slab.cell.lengths()
        slab_long_short_ratio=max(slab_length)/min(slab_length)
        if slab_long_short_ratio > 15:  
            calc=GPAW(xc=xc,
                h=h,
                symmetry = {'point_group': False},
                kpts=kpts,
                eigensolver=Davidson(3),
                mixer=Mixer(0.02, 5, 100),
                maxiter=500,
                occupations={'name': 'fermi-dirac','width': sw},
                poissonsolver={'dipolelayer': 'xy'})
        else:
            calc=GPAW(xc=xc,
                h=h,
                symmetry = {'point_group': False},
                kpts=kpts,
                eigensolver=Davidson(3),
                occupations={'name': 'fermi-dirac','width': sw},
                poissonsolver={'dipolelayer': 'xy'})            
        slab.set_calculator(calc)
        location=element+'/'+'surf'+'/'+struc+'/'+str(actual_layer)+'x1x1'
        opt.surf_relax(slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
        db_layer.write(slab,sim_layer=sim_layer,act_layer=actual_layer) #sim layer is different from the actual layers
        if iters>=2:
            fst=db_layer.get_atoms(id=iters-1)
            snd=db_layer.get_atoms(id=iters)
            trd=db_layer.get_atoms(id=iters+1)
            diff_primary=max(surf_e_calc(fst,snd,opt_bulk.get_potential_energy(),len(opt_bulk.get_tags())),surf_e_calc(fst,trd,opt_bulk.get_potential_energy(),len(opt_bulk.get_tags())))
            diff_second=surf_e_calc(snd,trd,opt_bulk.get_potential_energy(),len(opt_bulk.get_tags()))
            if temp_print==True:
                temp_output_printer(db_layer,iters,'act_layer',opt_bulk.get_potential_energy(),len(opt_bulk.get_tags()),rep_location)
        act_layer_ls.append(actual_layer)
        sim_layer_ls.append(sim_layer)
        iters+=1
        init_layer+=2 #change to one because the unit cell will generate 2 surfaces per layer
    
    
    if iters>=5:
        if diff_primary>rela_tol or diff_second>rela_tol:
            # Fiorentini and Methfessel relation (linear fit)
            with paropen(rep_location,'a') as f:
                parprint('Regular surface convergence failed.',file=f)
                parprint('Entering Fiorentini and Methfessel relation (linear fit) convergence test.',file=f)
            energy_slabs=[db_layer.get_atoms(i+1).get_potential_energy() for i in len(db_layer)]
            num_atoms=[db_layer.get(i+1).natoms for i in len(db_layer)]
            energy_bulk_fit=np.round(np.polyfit(num_atoms,energy_slabs,1)[0],decimals=5)
            fit_iters=2
            while (diff_primary>rela_tol or diff_second>rela_tol) and fit_iters <= 5:
                fst=db_layer.get_atoms(id=fit_iters-1)
                snd=db_layer.get_atoms(id=fit_iters)
                trd=db_layer.get_atoms(id=fit_iters+1)
                diff_primary=max(surf_e_calc(fst,snd,energy_bulk_fit,1),surf_e_calc(fst,trd,energy_bulk_fit,2))
                diff_second=surf_e_calc(snd,trd,energy_bulk_fit,1)
                if temp_print==True:
                    temp_output_printer(db_layer,iters,'act_layer',energy_bulk_fit,1,rep_location)
                fit_iters+=1
            if diff_primary>rela_tol or diff_second>rela_tol:
                with paropen(rep_location,'a') as f:
                    parprint("WARNING: Max Surface iterations reached! System may not be converged.",file=f)
                    parprint("Computation Suspended!",file=f)
                f.close()
                sys.exit()
            act_layer=act_layer_ls[fit_iters-3]
            sim_layer=sim_layer_ls[fit_iters-3]
            final_slab=db_layer.get_atoms(fit_iters-2) 
    act_layer=act_layer_ls[-3]
    sim_layer=sim_layer_ls[-3]
    final_slab=db_layer.get_atoms(len(db_layer)-2) 
    vac=final_slab.cell.lengths()[-1]-final_slab.positions[-1,2]
    db_final=connect('final_database'+'/'+'surf.db')
    id=db_final.reserve(name=element+'('+struc+')')
    if id is None:
        id=db_final.get(name=element+'('+struc+')').id
        db_final.update(id=id,atoms=final_slab,h=h,k_density=k_density,sw=sw,name=element+'('+struc+')',xc=xc,act_layer=act_layer,sim_layer=sim_layer,vac=vac)
    else:
        db_final.write(final_slab,id=id,name=element+'('+struc+')',h=h,k_density=k_density,sw=sw,xc=xc,act_layer=act_layer,sim_layer=sim_layer,vac=vac)
    with paropen(rep_location,'a') as f:
        parprint('Final Parameters:',file=f)
        parprint('\t'+'Simulated Layer: '+str(sim_layer),file=f)
        parprint('\t'+'Actual Layer: '+str(act_layer),file=f)
        parprint('\t'+'Vacuum length: '+str(vac)+'Ang',file=f)
        parprint('\t'+'Fixed layer: '+str(fix_layer),file=f)
        parprint('\t'+'xc: '+xc,file=f)
        parprint('\t'+'h: '+str(h),file=f)
        parprint('\t'+'k_density: '+str(k_density),file=f)
        parprint('\t'+'sw: '+str(sw),file=f)
    f.close()

def surf_e_calc(pre,post,bulk_e,bulk_num):
    #bulk_num=len(bulk.get_tags())
    #bulk_pot_e=bulk.get_potential_energy()
    opt_bulk_e=bulk_e/bulk_num
    pre_area=2*(pre.cell[0][0]*pre.cell[1][1])
    post_area=2*(post.cell[0][0]*post.cell[1][1])
    pre_e=pre.get_potential_energy()
    post_e=post.get_potential_energy()
    pre_num=float(re.findall(r'\d+',str(pre.symbols))[0])
    post_num=float(re.findall(r'\d+',str(post.symbols))[0])
    pre_surf_e=(1/pre_area)*(pre_e-pre_num*opt_bulk_e)
    post_surf_e=(1/post_area)*(post_e-post_num*opt_bulk_e)
    diff_surf_e=100*(abs((post_surf_e-pre_surf_e)/pre_surf_e))
    return diff_surf_e

def temp_output_printer(db,iters,key,bulk_e,bulk_num,location):
    fst_r=db.get(iters-1)
    snd_r=db.get(iters)
    trd_r=db.get(iters+1)
    with paropen(location,'a') as f:
        parprint('Optimizing parameter: '+key,file=f)
        parprint('\t'+'1st: '+str(fst_r[key])+' 2nd: '+str(snd_r[key])+' 3rd: '+str(trd_r[key])+'\n',file=f)
        parprint('\t'+'(2nd-1st)/1st: '+str(np.round(surf_e_calc(db.get_atoms(iters),db.get_atoms(iters-1),bulk_e,bulk_num),decimals=5))+'%',file=f)
        parprint('\t'+'(3nd-1st)/1st: '+str(np.round(surf_e_calc(db.get_atoms(iters+1),db.get_atoms(iters-1),bulk_e,bulk_num),decimals=5))+'%',file=f)
        parprint('\t'+'(3nd-2st)/2nd: '+str(np.round(surf_e_calc(db.get_atoms(iters+1),db.get_atoms(iters),bulk_e,bulk_num),decimals=5))+'%',file=f)
    f.close()