def make_mol_input_files(xyz_file,old_xyz_name,sample_file_name,total_charge=0,spin=0):
    inp_file_handle=open(sample_file_name,'r')
    xyz_name=xyz_file.split('/')[-1]
    element=xyz_name.split('.')[0]
    cid='_'.join(element.split('_')[-2:])
    output_file_handle=open('results'+'/'+cid+'/'+'mol.inp','w')
    for line in inp_file_handle:
        if old_xyz_name not in line:
            output_file_handle.write(line+'\n')
            continue
        if line.startswith('* xyzfile'):
            newline='* xyzfile'+' '+str(total_charge)+' '+str(2*spin+1)+' '
            output_file_handle.write(newline+xyz_file+'\n')
            continue
    output_file_handle.close()
    inp_file_handle.close()

