#!/usr/bin/env python

import pickle,subprocess,os,sys
import numpy as np
import pickle

# get work function with explicit QE
def get_wf_explicit(path):
    print(os.path.isfile(path+'/wf.out'))
    if os.path.isfile(path+'/wf.out'):
        f = open(path+'/wf.out','r')
        lines = f.readlines()
        f.close()
        return float(lines[0].rstrip())
    elif os.path.isfile(path+'/out.WF'):
        f = open(path+'/out.WF')
        lines = f.readlines()[0]
        return float(lines.split(',')[0].split('[')[-1])
    else:
        if os.path.isfile(path+'pot.xsf'):
            pot = pickle.load(open(path+'/wf_calc/pot.xsf','rb'),encoding='bytes')
            # load pickle file and take only the potential data
            pot = pot[b'data']

            # get planar average by taking mean across x and y axes
            pavg = np.mean(np.mean(pot,axis=0),axis=0)
        elif os.path.isfile(path+'potential.pickle'):
            potential = pickle.load(open(path+'potential.pickle', 'rb'), encoding='latin1')
            origin=potential[0]
            cell = potential[1]
            cell_x = cell[0][0] ; cell_y = cell[1][1] ; cell_z = cell[2][2]
            v = potential[2]
            n_x = len(v) # Returns the number of entries in the first subset (number of X entries)
            n_y = len(v[0]) # Returns the number of entries in the second subset (number of Y entries - choice of '0' was arbitrary)
            n_z = len(v[0][0]) # Returns the number of entries in the third subset (number of Z entries - choice of '0' was arbitrary)

            pavg = [] # Initializes a list for average potentials

            # Iterate over all space to average and process potential information
            for z in range(n_z):
                    temp = 0 # Will serve as a partial sum for averaging purposes
                    for x in range(n_x):
                            for y in range(n_y):
                                    temp+=v[x][y][z] # Sums the potential at all points in an X-Y plane
                    temp = temp / (n_x*n_y) # Divides by the total number of points
                    pavg.append( temp ) # Stores the value for average potential at this Z position

            # Translates the indices of points along the Z axis into real-space positions
            z_coordinate = np.linspace( 0 , cell_z , n_z )

        # assumes the system is positioned with (roughly) atoms in the center of the cell
        # with vacuum above and below
        lower_bound = int(9*len(pavg)/10)
        upper_bound = int(len(pavg))
        vac = np.mean(pavg[lower_bound:upper_bound])

        # grep log file for fermi level
        # using absolutely ridiculous python 3 subprocess command
        out = subprocess.Popen("grep 'the Fermi energy is' "+path+"/log | tail -n 1",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[-2]
        fermi = float(out.split()[-2])

        f = open(path+'/wf.out','w')
        f.writelines('%6f\n'%(vac-fermi))
        f.close()
        return vac-fermi

# get work function with Environ
def get_wf_implicit(path):
    # grep log file for fermi level
    # using absolutely ridiculous python 3 subprocess command
    out = subprocess.Popen("grep 'the Fermi energy is' "+path+"/log | tail -n 1",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[-2]
    fermi = float(out.split()[-2])

    try: 
        out = subprocess.Popen("grep ' due to the Gaussian-smeared nuclei' "+path+"/log | tail -n 1",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[-2]
    except IndexError:
        out = subprocess.Popen("grep ' due to the parabolic pbc-correction is' "+path+"/log | tail -n 1",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[-2]

    shift = float(out.split()[-2])

    return -1*(fermi+shift)

def get_wf_initial_implicit(path):
    # grep log file for fermi level
    # using absolutely ridiculous python 3 subprocess command
    out = subprocess.Popen("grep 'the Fermi energy is' "+path+"/log | tail -n 1",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[0]
    fermi = float(out.split()[-2])

    try: 
        out = subprocess.Popen("grep ' due to the Gaussian-smeared nuclei' "+path+"/log | tail -n 1",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[0]
    except IndexError:
        out = subprocess.Popen("grep ' due to the parabolic pbc-correction is' "+path+"/log | tail -n 1",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[1]

    shift = float(out.split()[-2])

    return -1*(fermi+shift)

def get_wf(path):
    # check if environ was used
    test = subprocess.Popen("grep 'Environ Module' "+path+"/log",shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    if  test == '':
        return get_wf_explicit(path)
    else:
        return get_wf_implicit(path)

#if len(sys.argv) > 1:
#    print(sys.argv)
#    print(get_wf(sys.argv[1]))
