from gpaw import GPAW,Mixer,Davidson
import bulk_autoconv as bulk_ac
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
## Set up the initial calculator
element='Cu_mp-30'
element_atom=bulk_ac.bulk_builder(element)
kpts=kdens2mp(element_atom) #default k_dens=3.5, even=True
calc=GPAW(xc='PBE',
            h=0.16,
            kpts=kpts,
            spinpol=False,
            maxiter=333,
            mixer=Mixer(0.05,5,50),
            eigensolver=Davidson(3),
            occupations={'name':'fermi-dirac','width':0.1})
bulk_ac.bulk_auto_conv(element,calc,rela_tol=10*10**(-3),temp_print=True)