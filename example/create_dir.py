#create directory before runing convergence test
import utils as ut
element_ls=['Cu','LiZn']
for element in element_ls:
    ut.create_dir(element,options=['bulk','surf'],
                surf_struc=['100','110','111'],
                optimized_parameters=['h','k','sw'],
                starting_layer=5,
                init_vac=8)
