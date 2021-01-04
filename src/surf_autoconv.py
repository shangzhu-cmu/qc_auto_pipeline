from gpaw import GPAW
from ase.constraints import FixAtoms
from ase.build import surface
from ase.io import write
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
def surf_auto_conv(element,struc,init_layer=5,vac=5,fix_layer=2,rela_tol=10*10**(-3),temp_print=True):
    # dispatcher = {'fcc100':build.fcc100,'fcc110':build.fcc110,'fcc111':build.fcc111,
    #               'bcc100':build.bcc100,'bcc110':build.bcc110,'bcc111':build.bcc111,
    #               }
    #m_ind=re.findall(r'\d+',struc)[0]
    m_ind=tuple(map(int,struc))
    rep_location=(element+'/'+'surf'+'/'+struc+'_results_report.txt')
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    try:
        db_bulk=connect('final_database'+'/'+'bulk.db')
    except:
        with paropen(rep_location,'a') as f:
            parprint('ERROR: No Optimized Bulk Object Found!',file=f)
            parprint('Surface Convergence Computation Suspended!',file=f)
        f.close()
        sys.exit()
    opt_bulk=db_bulk.get_atoms(name=element)
    xc=db_bulk.get(name=element).xc
    h=db_bulk.get(name=element).h
    k_density=db_bulk.get(name=element).k_density
    kpts=db_bulk.get(name=element).kpts
    sw=db_bulk.get(name=element).sw
    db_layer=connect(element+'/'+'surf'+'/'+'layer_converge.db')
    with paropen(rep_location,'a') as f:
        parprint('Initial Parameters:',file=f)
        parprint('\t'+'Materials: '+element,file=f)
        parprint('\t'+'xc: '+xc,file=f)
        parprint('\t'+'h: '+str(h),file=f)
        parprint('\t'+'k_density: '+str(k_density),file=f)
        parprint('\t'+'kpts: '+str(kpts),file=f)
        parprint('\t'+'sw: '+str(sw),file=f)
        parprint('\t'+'Miller Index: '+str(m_ind),file=f)
        parprint('\t'+'Layer: '+str(init_layer),file=f)
        parprint('\t'+'Vacuume length: '+str(vac)+'Ang',file=f)
        parprint('\t'+'Fixed layer: '+str(fix_layer),file=f)
        parprint('\t'+'rela_tol: '+str(rela_tol)+'eV',file=f)
    f.close()
    #first optimize the layers
    diff_primary=100
    diff_second=100
    iters=0
    layer_ls=[]
    area_rela_tol=0
    while (diff_primary>area_rela_tol or diff_second>area_rela_tol) and iters <= 6:
        slab = surface(opt_bulk, m_ind, layers=init_layer, vacuum=vac)
        fix_mask=slab.positions[:,2] <= np.unique(slab.positions[:,2])[fix_layer-1]
        slab.set_constraint(FixAtoms(mask=fix_mask))
        slab.set_pbc([1,1,0])
        kpts=kdens2mp(slab,kptdensity=k_density,even=True)
        calc=GPAW(xc=xc,
            h=h,
            symmetry = {'point_group': False},
            kpts=kpts,
            occupations={'name': 'fermi-dirac','width': sw},
            poissonsolver={'dipolelayer': 'xy'})
        slab.set_calculator(calc)
        location=element+'/'+'surf'+'/'+struc+'/'+str(init_layer)+'x1x1'
        opt.surf_relax(slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
        db_layer.write(slab,layers=init_layer)
        if iters>=2:
            area_rela_tol=cp.deepcopy(rela_tol)
            fst=db_layer.get_atoms(id=iters-1)
            snd=db_layer.get_atoms(id=iters)
            trd=db_layer.get_atoms(id=iters+1)
            diff_primary=max(surf_e_calc(fst,snd,opt_bulk),surf_e_calc(fst,trd,opt_bulk))
            diff_second=surf_e_calc(snd,trd,opt_bulk)
            if temp_print==True:
                temp_output_printer(db_layer,iters,'layers',opt_bulk,rela_tol,rep_location)
            area_average=np.mean([2*(fst.cell[0][0]*fst.cell[1][1]),
                                2*(snd.cell[0][0]*snd.cell[1][1]),
                                2*(trd.cell[0][0]*trd.cell[1][1])])
            area_rela_tol=area_rela_tol/area_average
        layer_ls.append(init_layer)
        iters+=1
        init_layer+=2
    if iters>=6:
        if diff_primary>area_rela_tol or diff_second>area_rela_tol:
            with paropen(rep_location,'a') as f:
                parprint("WARNING: Max LAYER iterations reached! System may not be converged.",file=f)
                parprint("Computation Suspended!",file=f)
            f.close()
            sys.exit()
    layer=layer_ls[-3]
    #second optimize the vaccum layer
    db_vac=connect(element+'/'+'surf'+'/'+'vac_converge.db')
    diff_primary=100
    diff_second=100
    iters=1
    vac_ls=[vac]
    db_vac.write(db_layer.get_atoms(len(db_layer)-2),vac=vac)
    area_rela_tol=0
    while (diff_primary>area_rela_tol or diff_second>area_rela_tol) and iters <= 5:
        slab = surface(opt_bulk, m_ind, layers=layer, vacuum=vac)
        fix_mask=slab.positions[:,2] <= np.unique(slab.positions[:,2])[fix_layer-1]
        slab.set_constraint(FixAtoms(mask=fix_mask))
        slab.set_pbc([1,1,0])
        kpts=kdens2mp(slab,kptdensity=k_density,even=True)
        calc=GPAW(xc=xc,
            h=h,
            symmetry = {'point_group': False},
            kpts=kpts,
            occupations={'name': 'fermi-dirac','width': sw},
            poissonsolver={'dipolelayer': 'xy'})
        slab.set_calculator(calc)
        location=element+'/'+'surf'+'/'+struc+'/'+'layer_optimized'+'/'+'vacuum_'+str(vac)
        opt.surf_relax(slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
        db_vac.write(slab,vac=vac)
        if iters>=2:
            area_rela_tol=cp.deepcopy(rela_tol)
            fst=db_vac.get_atoms(id=iters-1)
            snd=db_vac.get_atoms(id=iters)
            trd=db_vac.get_atoms(id=iters+1)
            diff_primary=max(surf_e_calc(fst,snd,opt_bulk),surf_e_calc(fst,trd,opt_bulk))
            diff_second=surf_e_calc(snd,trd,opt_bulk)
            if temp_print==True:
                temp_output_printer(db_layer,iters,'vac',opt_bulk,rela_tol,rep_location)
            area_average=np.mean([2*(fst.cell[0][0]*fst.cell[1][1]),
                                2*(snd.cell[0][0]*snd.cell[1][1]),
                                2*(trd.cell[0][0]*trd.cell[1][1])])
            area_rela_tol=area_rela_tol/area_average
        vac_ls.append(vac)
        iters+=1
        vac=int(vac+1)
    if iters>=5:
        if diff_primary>area_rela_tol or diff_second>area_rela_tol:
            with paropen(rep_location,'a') as f:
                parprint("WARNING: Max VACUUME iterations reached! System may not be converged.",file=f)
                parprint("Computation Suspended!",file=f)
            f.close()
            sys.exit()
    final_slab=db_vac.get_atoms(len(db_vac)-1)
    vac=vac_ls[-2]   
    db_final=connect('final_database'+'/'+'surf.db')
    id=db_final.reserve(name=element+struc)
    if id is None:
        id=db_final.get(name=element+struc).id
        db_final.update(id=id,atoms=final_slab,h=h,k_density=k_density,sw=sw,name=element+struc,xc=xc,layer=layer,vac=vac)
    else:
        db_final.write(final_slab,id=id,name=element+struc,h=h,k_density=k_density,sw=sw,xc=xc,layer=layer,vac=vac)
    with paropen(rep_location,'a') as f:
        parprint('Final Parameters:',file=f)
        parprint('\t'+'xc: '+xc,file=f)
        parprint('\t'+'h: '+str(h),file=f)
        parprint('\t'+'k_density: '+str(k_density),file=f)
        parprint('\t'+'sw: '+str(sw),file=f)
        parprint('\t'+'Layer: '+str(init_layer),file=f)
        parprint('\t'+'Vacuume length: '+str(vac)+'Ang',file=f)
        parprint('\t'+'Fixed layer: '+str(fix_layer),file=f)
        parprint(time_ls,file=f)
        parprint(vac_ls,file=f)
    f.close()

def surf_e_calc(pre,post,bulk):
    bulk_num=len(bulk.get_tags())
    bulk_pot_e=bulk.get_potential_energy()
    opt_bulk_e=bulk_pot_e/bulk_num
    pre_area=2*(pre.cell[0][0]*pre.cell[1][1])
    post_area=2*(post.cell[0][0]*post.cell[1][1])
    pre_e=pre.get_potential_energy()
    post_e=post.get_potential_energy()
    pre_num=float(re.findall(r'\d+',str(pre.symbols))[0])
    post_num=float(re.findall(r'\d+',str(post.symbols))[0])
    pre_surf_e=(1/pre_area)*(pre_e-pre_num*opt_bulk_e)
    post_surf_e=(1/post_area)*(post_e-post_num*opt_bulk_e)
    diff_surf_e=abs(post_surf_e-pre_surf_e)
    return diff_surf_e

def temp_output_printer(db,iters,key,bulk,rela_tol,location):
    fst_r=db.get(iters-1)
    snd_r=db.get(iters)
    trd_r=db.get(iters+1)
    area_average=np.mean([2*(fst_r.cell[0][0]*fst_r.cell[1][1]),
                        2*(snd_r.cell[0][0]*snd_r.cell[1][1]),
                        2*(trd_r.cell[0][0]*trd_r.cell[1][1])])
    rela_tol/area_average
    with paropen(location,'a') as f:
        parprint('Optimizing parameter: '+key,file=f)
        parprint('Relative Tolerance: '+str(np.round(rela_tol/area_average,decimals=5))+'eV/Ang^2',file=f)
        parprint('\t'+'1st: '+str(fst_r[key])+' 2nd: '+str(snd_r[key])+' 3rd: '+str(trd_r[key])+'\n',file=f)
        parprint('\t'+'2nd-1st: '+str(np.round(surf_e_calc(db.get_atoms(iters),db.get_atoms(iters-1),bulk),decimals=5))+'eV/Ang^2',file=f)
        parprint('\t'+'3nd-1st: '+str(np.round(surf_e_calc(db.get_atoms(iters+1),db.get_atoms(iters-1),bulk),decimals=5))+'eV/Ang^2',file=f)
        parprint('\t'+'3nd-2st: '+str(np.round(surf_e_calc(db.get_atoms(iters+1),db.get_atoms(iters),bulk),decimals=5))+'eV/Ang^2',file=f)
    f.close()