import bulk_autoconv as bulk_ac
import surf_autoconv as surf_ac
element_ls=['Cu','LiZn'] #Cu can be built using ase.build.bulk; LiZn can be built using .cif
a_ls=[3.65,None]
struc_ls=['fcc',None]
h=0.16
k_dens=2.5
sw=0.1
rela_tol=15*10**(-3) #eV
cif_ls=[True, False]
option=True
for i in range(len(element_ls)):
    bulk_ac.bulk_auto_conv(element_ls[i],
                            a0=a_ls[i],
                            struc=struc_ls[i],
                            h=h,
                            k_density=k_dens,
                            xc='PBE',
                            sw=sw,
                            rela_tol=rela_tol,
                            cif=cif_ls[i],
                            temp_print=option)
for i in range(len(element_ls)):
    for m_ind in ['110','100','111']:
        surf_ac.surf_auto_conv(element_ls[i],
                                m_ind,
                                init_layer=5,
                                vac=8,
                                fix_layer=2,
                                rela_tol=rela_tol,
                                temp_print=option)


    

