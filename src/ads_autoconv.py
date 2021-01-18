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
import copy as cp
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
from autocat import adsorption
import glob
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

def ads_auto_conv(element,struc,ads,ads_pot_e,ads_site=['ontop','hollow','bridge'],fix_layer=2,rela_tol=5,temp_print=True): #changed the rela_tol to percentage 
    #convert str ind to tuple
    m_ind=tuple(map(int,struc))

    #set up the workspace
    code_dir=os.getcwd() #get the current working dir
    struc_dir=element+'/'+'ads'+'/'+struc #get the structure dir
    os.chdir(code_dir+'/'+struc_dir) #change the working dir to structure dir "./element/ads/struc"

    #create report
    rep_location=code_dir+'/'+element+'/'+'ads'+'/'+struc+'_results_report.txt' 
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    
    #connect to the surface database to get the parameters for calculation
    #connect to bulk database to get the optimized bulk for clean slab generation preparation
    db_surf=connect(code_dir+'/'+'final_database'+'/'+'surf.db')
    opt_bulk=connect(code_dir+'/'+'final_database'+'/'+'bulk.db').get_atoms(name=element)
    pymatgen_bulk=AseAtomsAdaptor.get_structure(opt_bulk)

    #get out the parameters for calculation from surfae database
    xc=db_surf.get(name=element+'('+struc+')').xc
    h=db_surf.get(name=element+'('+struc+')').h
    k_density=db_surf.get(name=element+'('+struc+')').k_density
    sw=db_surf.get(name=element+'('+struc+')').sw
    act_init_layer=db_surf.get(name=element+'('+struc+')').act_layer
    sim_init_layer=db_surf.get(name=element+'('+struc+')').sim_layer
    vac=db_surf.get(name=element+'('+struc+')').vac
    
    #print out parameters
    with paropen(rep_location,'a') as f:
        parprint('Initial Parameters:',file=f)
        parprint('\t'+'Materials: '+element,file=f)
        parprint('\t'+'Miller Index: '+str(m_ind),file=f)
        parprint('\t'+'Adsorbate: '+ads,file=f)
        parprint('\t'+'Simulated Layer: '+str(sim_init_layer),file=f)
        parprint('\t'+'Actual Layer: '+str(act_init_layer),file=f)
        parprint('\t'+'Vacuum Length: '+str(vac)+'Ang',file=f)
        parprint('\t'+'Fixed Layer: '+str(fix_layer),file=f)
        parprint('\t'+'xc: '+xc,file=f)
        parprint('\t'+'h: '+str(h),file=f)
        parprint('\t'+'k_density: '+str(k_density),file=f)
        parprint('\t'+'sw: '+str(sw),file=f)
        parprint('\t'+'rela_tol: '+str(rela_tol)+'%',file=f)
    f.close()

    #optimize the layers
    ##connect to adsorbed layer database
    ##connect to clean layer database
    db_ads_slab=connect(code_dir+'/'+element+'/'+'ads'+'/'+struc+'/'+'layer_converge.db')
    db_slab_clean=connect(code_dir+'/'+element+'/'+'surf'+'/'+struc+'/'+'layer_converge.db')
    
    ##setup parameters for the following convergence test
    avail_slab_num=len(db_slab_clean)#number of previous tests (clean slab)
    diff_primary=100 
    diff_second=100
    iters=0
    act_layer_ls=[]
    sim_layer_ls=[]
    ads_file_loc=code_dir+'/'+element+'/'+'ads'+'/'+struc#save the location of adsorption folder (./element/ads/struc)
    
    while (diff_primary>rela_tol or diff_second>rela_tol) and iters < avail_slab_num:
        #glob all the .traj file in the './element/nx1x1/Li/**/' folder
        act_init_layer=db_slab_clean.get(iters+1).act_layer
        sim_init_layer=db_slab_clean.get(iters+1).sim_layer
        fil=glob.glob(ads_file_loc+'/'+str(act_init_layer)+'x1x1'+'/'+'adsorbates/Li/**/**/*.traj',recursive=False)
        ads_dict={} #create dictionary for saving the adsorption energy
        for file_loc in fil: 
            ads_slab = read(file_loc)
            kpts=kdens2mp(ads_slab,kptdensity=k_density,even=True)
            calc=GPAW(xc='PBE',
                    h=h,
                    symmetry = {'point_group': False},
                    kpts=kpts,
                    occupations={'name': 'fermi-dirac','width': sw},
                    poissonsolver={'dipolelayer': 'xy'})
            ads_slab.set_calculator(calc)
            location='/'.join(file_loc.split('/')[:-1])
            opt.surf_relax(ads_slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
            ads_dict[location]=ads_slab.get_potential_energy()-(db_slab_clean.get_atoms(iters+1).get_potential_energy()+ads_pot_e)
        ads_dict_sorted=sorted(ads_dict,key=ads_dict.get)
        lowest_ads_e_slab=read(ads_dict_sorted[0])
        db_ads_slab.write(lowest_ads_e_slab,act_layer=act_init_layer,sim_layer=sim_init_layer)
        #enter the convergence test sequence
        if iters>=2:
            fst_ads=db_ads_slab.get_atoms(id=iters-1)
            snd_ads=db_ads_slab.get_atoms(id=iters)
            trd_ads=db_ads_slab.get_atoms(id=iters+1)
            fst_clean=db_slab_clean.get_atoms(id=iters-1)
            snd_clean=db_slab_clean.get_atoms(id=iters)
            trd_clean=db_slab_clean.get_atoms(id=iters+1)
            ads_e_perct_diff_12=100*abs((ads_e_calc(fst_ads,fst_clean,ads_pot_e)-ads_e_calc(snd_ads,snd_clean,ads_pot_e))/ads_e_calc(fst_ads,fst_clean,ads_pot_e))
            ads_e_perct_diff_13=100*abs((ads_e_calc(fst_ads,fst_clean,ads_pot_e)-ads_e_calc(trd_ads,trd_clean,ads_pot_e))/ads_e_calc(fst_ads,fst_clean,ads_pot_e))
            ads_e_perct_diff_23=100*abs((ads_e_calc(snd_ads,snd_clean,ads_pot_e)-ads_e_calc(trd_ads,trd_clean,ads_pot_e))/ads_e_calc(snd_ads,snd_clean,ads_pot_e))
            diff_primary=max(ads_e_perct_diff_12,ads_e_perct_diff_13)
            diff_second=ads_e_perct_diff_23
            if temp_print==True:
                temp_output_printer(db_ads_slab,db_slab_clean,iters,'act_layer',ads_pot_e,rep_location)
        act_layer_ls.append(act_init_layer)
        sim_layer_ls.append(sim_init_layer)
        iters+=1
    
    #entering the backup running sequence
    if (diff_primary>rela_tol or diff_second>rela_tol):
        #write warning message
        with paropen(rep_location,'a') as f:
            parprint('Adsorption energy did not converge with the pre-computed slabs.',file=f)
            parprint('WARNING: New computation from clean slab can take a long time.',file=f)
            parprint('Consider terminate the computation if the relative percentage differences are around the desire range', file=f)
        #entering while loop the iterations will stop when it hits 6
        sim_init_layer+=1
        act_init_layer+=2
        # clean_slab = surface(opt_bulk, m_ind, layers=sim_init_layer, vacuum=vac)
        slabgen = SlabGenerator(pymatgen_bulk, m_ind, sim_init_layer, sim_init_layer*2, center_slab=True, lll_reduce=True, in_unit_planes=True)
        slabs=slabgen.get_slabs() 
        slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
        clean_slab=AseAtomsAdaptor.get_atoms(slabs_symmetric[0]) #convert to ase structure
        # orthogonality=clean_slab.get_cell_lengths_and_angles()[3:5]==90
        # if not np.all(orthogonality):
        #     clean_slab=surface(opt_bulk, m_ind, layers=sim_init_layer, vacuum=vac)
        #     orthogonality=clean_slab.get_cell_lengths_and_angles()[3:5]==90    
        #     if not np.all(orthogonality):
        #         with paropen(rep_location,'a') as f:
        #             parprint('ERROR: Cannot create surface with orthogonality.',file=f)
        #             parprint('Computation Suspended!',file=f)
        #             sys.exit()
        actual_layer=len(np.unique(np.round(clean_slab.positions[:,2],decimals=4)))
        while (diff_primary>rela_tol or diff_second>rela_tol) and iters < 6:
            #clean slab optimization
            while actual_layer != act_init_layer:
                sim_init_layer+=1
                #clean_slab=surface(opt_bulk, m_ind, layers=sim_init_layer,vacuum=vac)
                slabgen = SlabGenerator(pymatgen_bulk, m_ind, sim_init_layer, sim_init_layer*2, center_slab=True, lll_reduce=True, in_unit_planes=True)
                slabs=slabgen.get_slabs() 
                slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
                clean_slab=AseAtomsAdaptor.get_atoms(slabs_symmetric[0]) #convert to ase structure
                # orthogonality=clean_slab.get_cell_lengths_and_angles()[3:5]==90
                # if not np.all(orthogonality):
                #     clean_slab=surface(opt_bulk, m_ind, layers=sim_init_layer, vacuum=vac)
                #     orthogonality=clean_slab.get_cell_lengths_and_angles()[3:5]==90    
                #     if not np.all(orthogonality):
                #         with paropen(rep_location,'a') as f:
                #             parprint('ERROR: Cannot create surface with orthogonality.',file=f)
                #             parprint('Computation Suspended!',file=f)
                #             sys.exit()
                actual_layer=len(np.unique(np.round(clean_slab.positions[:,2],decimals=4)))
                if actual_layer > act_init_layer:
                    with paropen(rep_location,'a') as f:
                        parprint('ERROR: Actual number of layers is greater than the desired number of layers.',file=f)
                        parprint('\t'+'Actual Layer: '+str(actual_layer),file=f)
                        parprint('\t'+'Desired Layer: '+str(act_init_layer),file=f)
                        parprint('Computation Suspended!',file=f)
                        sys.exit() 
            current_vac=clean_slab.cell.lengths()[-1]-clean_slab.positions[-1,2]
            while current_vac < vac:
                clean_slab.center(vacuum=vac,axis=2)
                # vac_init_layer+=1
                # slabgen = SlabGenerator(pymatgen_bulk, m_ind, sim_init_layer, vac_init_layer, center_slab=True, lll_reduce=True, in_unit_planes=True)
                # slabs=slabgen.get_slabs() 
                # slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
                # clean_slab=AseAtomsAdaptor.get_atoms(slabs_symmetric[0]) #convert to ase structure
                current_vac=clean_slab.cell.lengths()[-1]-clean_slab.positions[-1,2]                    
            fix_mask=np.round(clean_slab.positions[:,2],decimals=4) <= np.unique(np.round(clean_slab.positions[:,2],decimals=4))[fix_layer-1]
            clean_slab.set_constraint(FixAtoms(mask=fix_mask))
            clean_slab.set_pbc([1,1,0])
            kpts=kdens2mp(clean_slab,kptdensity=k_density,even=True)
            calc=GPAW(xc=xc,
                    h=h,
                    symmetry = {'point_group': False},
                    kpts=kpts,
                    occupations={'name': 'fermi-dirac','width': sw},
                    poissonsolver={'dipolelayer': 'xy'})
            clean_slab.set_calculator(calc)
            location=code_dir+'/'+element+'/'+'surf'+'/'+struc+'/'+str(actual_layer)+'x1x1'
            opt.surf_relax(clean_slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
            db_slab_clean.write(clean_slab,sim_layer=sim_init_layer,act_layer=actual_layer)
            #create the adsorption site file using autocat 
            ads_creater(clean_slab,ads,actual_layer,ads_site,ads_file_loc)
            os.chdir(code_dir+'/'+struc_dir)
            fil=glob.glob(ads_file_loc+'/'+str(actual_layer)+'x1x1'+'/'+'adsorbates/Li/**/**/*.traj',recursive=False)
            ads_dict={}
            for file_loc in fil: 
                ads_slab = read(file_loc)
                kpts=kdens2mp(ads_slab,kptdensity=k_density,even=True)
                calc=GPAW(xc='PBE',
                        h=h,
                        symmetry = {'point_group': False},
                        kpts=kpts,
                        occupations={'name': 'fermi-dirac','width': sw},
                        poissonsolver={'dipolelayer': 'xy'})
                ads_slab.set_calculator(calc)
                location='/'.join(file_loc.split('/')[:-1])
                opt.surf_relax(ads_slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
                ads_dict[location]=ads_slab.get_potential_energy()-(db_slab_clean.get_atoms(iters+1).get_potential_energy()+ads_pot_e)
            ads_dict_sorted=sorted(ads_dict,key=ads_dict.get)
            lowest_ads_e_slab=read(ads_dict_sorted[0])
            db_ads_slab.write(lowest_ads_e_slab,act_layer=actual_layer,sim_layer=sim_init_layer)
            #enter the convergence test sequence
            fst_ads=db_ads_slab.get_atoms(id=iters-1)
            snd_ads=db_ads_slab.get_atoms(id=iters)
            trd_ads=db_ads_slab.get_atoms(id=iters+1)
            fst_clean=db_slab_clean.get_atoms(id=iters-1)
            snd_clean=db_slab_clean.get_atoms(id=iters)
            trd_clean=db_slab_clean.get_atoms(id=iters+1)
            ads_e_perct_diff_12=100*abs((ads_e_calc(fst_ads,fst_clean,ads_pot_e)-ads_e_calc(snd_ads,snd_clean,ads_pot_e))/ads_e_calc(fst_ads,fst_clean,ads_pot_e))
            ads_e_perct_diff_13=100*abs((ads_e_calc(fst_ads,fst_clean,ads_pot_e)-ads_e_calc(trd_ads,trd_clean,ads_pot_e))/ads_e_calc(fst_ads,fst_clean,ads_pot_e))
            ads_e_perct_diff_23=100*abs((ads_e_calc(snd_ads,snd_clean,ads_pot_e)-ads_e_calc(trd_ads,trd_clean,ads_pot_e))/ads_e_calc(snd_ads,snd_clean,ads_pot_e))
            diff_primary=max(ads_e_perct_diff_12,ads_e_perct_diff_13)
            diff_second=ads_e_perct_diff_23
            if temp_print==True:
                temp_output_printer(db_ads_slab,db_slab_clean,iters,'act_layer',ads_pot_e,rep_location)
            act_layer_ls.append(actual_layer)
            sim_layer_ls.append(sim_init_layer)
            act_init_layer+=2
            iters+=1

    if iters>=5:
        if diff_primary>rela_tol or diff_second>rela_tol:
            with paropen(rep_location,'a') as f:
                parprint("WARNING: Max Adsorption iterations reached! System may not be converged.",file=f)
                parprint("Computation Suspended!",file=f)
            f.close()
            sys.exit()

    final_slab_ads=db_ads_slab.get_atoms(len(db_ads_slab)-2)
    final_slab=db_slab_clean.get_atoms(len(db_slab_clean)-2)
    act_layer=act_layer_ls[-3]
    sim_layer=sim_layer_ls[-3]
    db_final_full=connect(code_dir+'/'+'final_database'+'/'+'full_ads.db')
    db_final_clean=connect(code_dir+'/'+'final_database'+'/'+'clean_ads.db')
    id=db_final_full.reserve(name=element+'('+struc+')')
    if id is None:
        id=db_final_full.get(name=element+'('+struc+')').id
        db_final_full.update(id=id,atoms=final_slab_ads,h=h,k_density=k_density,sw=sw,name=element+'('+struc+')',xc=xc,act_layer=act_layer,sim_layer=sim_layer,vac=vac)
    else:
        db_final_full.write(final_slab_ads,id=id,h=h,k_density=k_density,sw=sw,name=element+'('+struc+')',xc=xc,act_layer=act_layer,sim_layer=sim_layer,vac=vac)
    id=db_final_clean.reserve(name=element+'('+struc+')')
    if id is None:
        id=db_final_clean.get(name=element+'('+struc+')').id
        db_final_clean.update(id=id,atoms=final_slab,h=h,k_density=k_density,sw=sw,name=element+'('+struc+')',xc=xc,act_layer=act_layer,sim_layer=sim_layer,vac=vac)
    else:
        db_final_clean.write(final_slab,id=id,h=h,k_density=k_density,sw=sw,name=element+'('+struc+')',xc=xc,act_layer=act_layer,sim_layer=sim_layer,vac=vac)
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

def ads_creater(slab,atom,act_init_layer,site,file_loc):
    dir_loc=file_loc+'/'+str(act_init_layer)+'x1x1'
    if world.rank==0:
        os.makedirs(dir_loc,exist_ok=True)
    os.chdir(dir_loc)
    adsorption.generate_rxn_structures(slab,ads=atom,site_type=site,write_to_disk=True)

    

def ads_e_calc(full,clean,ads_pot_e):
    ads_e=full.get_potential_energy()-(clean.get_potential_energy()+ads_pot_e)
    return ads_e

def temp_output_printer(db_layer,db_slab_clean,iters,key,ads_pot_e,location):
    fst_ads=db_layer.get_atoms(id=iters-1)
    snd_ads=db_layer.get_atoms(id=iters)
    trd_ads=db_layer.get_atoms(id=iters+1)
    fst_clean=db_slab_clean.get_atoms(id=iters-1)
    snd_clean=db_slab_clean.get_atoms(id=iters)
    trd_clean=db_slab_clean.get_atoms(id=iters+1)
    ads_e_perct_diff_12=100*abs((ads_e_calc(fst_ads,fst_clean,ads_pot_e)-ads_e_calc(snd_ads,snd_clean,ads_pot_e))/ads_e_calc(fst_ads,fst_clean,ads_pot_e))
    ads_e_perct_diff_13=100*abs((ads_e_calc(fst_ads,fst_clean,ads_pot_e)-ads_e_calc(trd_ads,trd_clean,ads_pot_e))/ads_e_calc(fst_ads,fst_clean,ads_pot_e))
    ads_e_perct_diff_23=100*abs((ads_e_calc(snd_ads,snd_clean,ads_pot_e)-ads_e_calc(trd_ads,trd_clean,ads_pot_e))/ads_e_calc(snd_ads,snd_clean,ads_pot_e))
    with paropen(location,'a') as f:
        parprint('Optimizing parameter: '+key,file=f)
        parprint('\t'+'1st: '+str(db_layer.get(iters-1)[key])+' 2nd: '+str(db_layer.get(iters)[key])+' 3rd: '+str(db_layer.get(iters+1)[key])+'\n',file=f)
        parprint('\t'+'(2nd-1st)/1st: '+str(np.round(ads_e_perct_diff_12,decimals=5))+'%',file=f)
        parprint('\t'+'(3nd-1st)/1st: '+str(np.round(ads_e_perct_diff_13,decimals=5))+'%',file=f)
        parprint('\t'+'(3nd-2st)/2nd: '+str(np.round(ads_e_perct_diff_23,decimals=5))+'%',file=f)
    f.close()