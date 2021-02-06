from gpaw import GPAW
import glob
import optimizer as opt
import numpy as np 
from ase.parallel import parprint,world,paropen
import os
from ase.db import connect
import sys
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
from ase.io import read,write
from gpaw import Davidson
from gpaw import Mixer,MixerDif
def ads_auto_select(element,
                struc,
                ads,
                ads_pot_e,
                maxiter=333,
                beta=0.05,
                nmaxold=5,
                weight=50.0,
                temp_print=True,
                size='1x1'):

    #convert str ind to tuple
    m_ind=tuple(map(int,struc))

    #set up the workspace
    code_dir=os.getcwd() #get the current working dir
    #struc_dir=element+'/'+'ads'+'/'+struc #get the structure dir
    # 

    #create report
    rep_location=code_dir+'/'+element+'/'+'ads'+'/'+struc+'_results_report.txt' 
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    if not os.path.isfile('final_database/surf.db'):
        with paropen(rep_location,'a') as f:
            parprint('ERROR: surf database has not been established!',file=f)
            parprint('Adsorption Computation Suspended!',file=f)
        f.close()
        sys.exit()
    else:
        db_surf=connect('final_database'+'/'+'surf.db')
        db_bulk=connect('final_database'+'/'+'bulk.db')
        try:
            opt_slab=db_surf.get_atoms(name=element+'('+struc+')')
        except:
            with paropen(rep_location,'a') as f:
                parprint('ERROR: No Optimized Surf Object Found!',file=f)
                parprint('Adsorption Computation Suspended!',file=f)
            f.close()
            sys.exit()
    #connect to the surface database to get the parameters for calculation

    xc=db_surf.get(name=element+'('+struc+')').xc
    h=db_surf.get(name=element+'('+struc+')').h
    k_density=db_surf.get(name=element+'('+struc+')').k_density
    sw=db_surf.get(name=element+'('+struc+')').sw
    spin=db_bulk.get(name=element).spin
    if spin:
        magmom=db_surf.get(name=element+'('+struc+')').magmom

    with paropen(rep_location,'a') as f:
        parprint('Initial Parameters:',file=f)
        parprint('\t'+'Materials: '+element,file=f)
        parprint('\t'+'Miller Index: '+str(m_ind),file=f)
        parprint('\t'+'Adsorbate: '+str(ads),file=f)
        parprint('\t'+'xc: '+xc,file=f)
        parprint('\t'+'h: '+str(h),file=f)
        parprint('\t'+'k_density: '+str(k_density),file=f)
        parprint('\t'+'sw: '+str(sw),file=f)
        parprint('\t'+'spin polarized: '+str(spin),file=f)
        if spin:
            parprint('\t'+'magmom: '+str(magmom),file=f)
    f.close()

    ads_db=connect('final_database/ads'+str(size)+'.db')
    
    ads_file_loc=code_dir+'/'+element+'/'+'ads'+'/'+struc
    fils=glob.glob(ads_file_loc+'/'+'adsorbates/Li/**/**/*.traj',recursive=False)
    ads_dict={}
    with paropen(rep_location,'a') as f:
        parprint('Ads Site(Ang)\t\t\tAds Energy(eV)',file=f)
    for file_loc in fils:
        ads_slab=read(file_loc)
        kpts=kdens2mp(ads_slab,kptdensity=k_density,even=True)
        slab_length=ads_slab.cell.lengths()
        slab_long_short_ratio=max(slab_length)/min(slab_length)
        if spin:
            slab_formula=ads_slab.get_chemical_symbols()
            magmom_ls=np.mean(magmom)*np.ones(len(slab_formula))
            magmom_ls[slab_formula.index(ads)]=0
            ads_slab.set_initial_magnetic_moments(magmom_ls)
        if slab_long_short_ratio > 15:  
            if spin:
                calc=GPAW(xc=xc,
                    h=h,
                    symmetry = {'point_group': False},
                    kpts=kpts,
                    eigensolver=Davidson(3),
                    mixer=MixerDif(np.round(beta/2,decimals=2), nmaxold, weight*2),
                    spinpol=spin,
                    maxiter=maxiter,
                    occupations={'name': 'fermi-dirac','width': sw},
                    poissonsolver={'dipolelayer': 'xy'})
            else:
                calc=GPAW(xc=xc,
                    h=h,
                    symmetry = {'point_group': False},
                    kpts=kpts,
                    eigensolver=Davidson(3),
                    mixer=Mixer(np.round(beta/2,decimals=2), nmaxold, weight*2),
                    maxiter=maxiter,
                    occupations={'name': 'fermi-dirac','width': sw},
                    poissonsolver={'dipolelayer': 'xy'})
        else:
            if spin:
                calc=GPAW(xc=xc,
                    h=h,
                    symmetry = {'point_group': False},
                    kpts=kpts,
                    spinpol=spin,
                    mixer=MixerDif(beta,nmaxold,weight),
                    maxiter=maxiter,
                    eigensolver=Davidson(3),
                    occupations={'name': 'fermi-dirac','width': sw},
                    poissonsolver={'dipolelayer': 'xy'})   
            else:
                calc=GPAW(xc=xc,
                    h=h,
                    symmetry = {'point_group': False},
                    kpts=kpts,
                    spinpol=spin,
                    mixer=Mixer(beta,nmaxold,weight),
                    maxiter=maxiter,
                    eigensolver=Davidson(3),
                    occupations={'name': 'fermi-dirac','width': sw},
                    poissonsolver={'dipolelayer': 'xy'})                            
        ads_slab.set_calculator(calc)
        location='/'.join(file_loc.split('/')[:-1])
        opt.surf_relax(ads_slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
        ads_dict[location]=ads_slab.get_potential_energy()-(opt_slab.get_potential_energy()+ads_pot_e)
        if temp_print:
            with paropen(rep_location,'a') as f:
                parprint(str(file_loc.split('/')[-2])+'\t\t\t'+str(np.round(ads_dict[location],decimals=5)),file=f)
    ads_dict_sorted=sorted(ads_dict,key=ads_dict.get)
    lowest_ads_e_slab=read(ads_dict_sorted[0]+'/slab.traj')
    ads_db.write(lowest_ads_e_slab)
    with paropen(rep_location,'a') as f:
        parprint('Computation Complete. Selected ads site is: '+ads_dict_sorted[0].split('/')[-1],file=f)