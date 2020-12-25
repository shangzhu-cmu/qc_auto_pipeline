import bulk_autoconv as bulk_ac
element='Cu'
bulk_ac.bulk_auto_conv(element,h=0.18, k=6,xc='PBE',sw=0.1, rela_tol=5*10**(-3))
