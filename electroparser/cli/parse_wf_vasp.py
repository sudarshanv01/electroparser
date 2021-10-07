#!/usr/bin/env python
import os,sys,subprocess
import electroparser.cli.NewPotentialModule as pot
from ase.io import read
import numpy as np

def get_wf_implicit(path, neb=False):
    out1 = subprocess.Popen('grep fermi '+path+'/OUTCAR | tail -n 1',shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    fermi = float(out1.split()[2])
    if neb == False:
        try:
            out2 = subprocess.Popen('grep FERMI_SHIFT '+path+'/vasp.out | tail -n 1',shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            shift = float(out2.split(' = ')[-1])
        except ValueError:
            out2 = subprocess.Popen('grep FERMI_SHIFT '+path+'/opt.log | tail -n 1',shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            shift = float(out2.split(' = ')[-1])
    elif neb == True:
        # if it a neb take the fermi shift from the directory above it
        try:
            out2 = subprocess.Popen('grep FERMI_SHIFT '+path+'../vasp.out | tail -n 1',shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            shift = float(out2.split(' = ')[-1])
        except ValueError:
            out2 = subprocess.Popen('grep FERMI_SHIFT '+path+'../slurm* | tail -n 1',shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            shift = float(out2.split(' = ')[-1])


    return -1*(fermi+shift)
def get_fermi_shift(path):
    out2 = subprocess.Popen('grep FERMI_SHIFT '+path+'/vasp.out | tail -n 1',shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    shift = float(out2.split(' = ')[-1])
    return shift
def get_wf_explicit(path):
    if os.path.isfile(path+'/wf.out'):
        f = open(path+'/wf.out')
        lines = f.readlines()
        return float(lines[0].rstrip())

    out = 'planar.dat'
    vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density(path+'/LOCPOT')
    vector_a,vector_b,vector_c,av,bv,cv = pot.matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ)
    #------------------------------------------------------------------
    ## POTENTIAL
    potential = pot.planar_average(grid_pot,NGX,NGY,NGZ)
    atoms = read(path+'/POSCAR')
    cell_z = atoms.cell[2][2]
    n_z = len(potential)
    z_coord = np.linspace(0,cell_z,n_z)

    # from ase, fix later
    # position_array = atoms.get_scaled_positions()[...,2]
    # position_array.sort()
    # diffs = np.diff(position_array)
    # diffs = np.append(diffs,position_array[0]+1-position_array[-1])
    # max_diff_index = np.argmax(diffs)
    # if max_diff_index == len(position_array)-1:
        # vacuum_pos = (position_array[0] +1 - position_array[-1])/2.0 % 1
    # else:
        # vacuum_pos = (position_array[max_diff_index] + position_array[max_diff_index + 1])/2.
    
    # # vacuum energy
    # vac_e = potential[np.abs(np.array(potential)[...,0]-vacuum_pos).argmin()][1]
    psi2_sum = 0
    range2 = range(int(8*n_z/10),int(9*n_z/10))
    for n in range2:
        psi2_sum += potential[n]
    psi2 = psi2_sum/len(range2)

    out1 = subprocess.Popen('grep fermi '+path+'/OUTCAR | tail -n 1',shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    fermi = float(out1.split()[2])
    psi2 -= fermi
    f = open(path + '/wf.out','w')
    f.writelines([str(psi2)+'\n'])
    f.close()
    return psi2

# this should be executed in the directory you want to get the WF.
# first, find out if it's an implicit solvent calculation
# print 'why is the function here'
# try:
    # out = subprocess.check_output('grep LSOL OUTCAR',shell=True).split()
    # if out[2] == 'T':
        # print get_wf_implicit('.')
    # else:
        # print get_wf_explicit()
# except subprocess.CalledProcessError as e:
    # if e.cmd == "grep LSOL OUTCAR":
        # print get_wf_explicit('.')

