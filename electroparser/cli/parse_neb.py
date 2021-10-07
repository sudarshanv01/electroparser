#!/usr/bin/python

from ase.io import read, write
from ase.db import connect
import numpy as np
import os, sys
from pathlib import Path
sys.path.append('/home/cat/vijays/tools/scripts')
import get_wf as qe_wf
import get_wf_vasp as vasp_wf
import re
from ase import units
import glob
import subprocess
from ase.vibrations import Vibrations
from glob import glob

###############################################
# Script to parse data from any set of folders and put it into 
# an ASE database
# To do:
# 1. Pick up the state from anywhere
# 2. Add in frequencies
# 3. Add in support for different problems
# 4. Get all environ keys
###############################################

def get_atoms(homedir, image):
    qe = Path(homedir + '/run_aseQE.py')
    vasp = Path(homedir + '/OUTCAR')
    
    if vasp.is_file():
        # gets the atoms object for each neb image
        atoms = read(homedir +  '/OUTCAR')
    elif qe.is_file():
        # Returns the neb_trajectory based on image number
    #    flist = glob(homedir + '/neb*.traj')
        atoms = read(homedir + '/neb' + str(image) +  '.traj')
    return atoms


        
def get_wf_parsed(basedir, images):
    # Returns the WF 
    # If there is none, retruns nan

    vasp = Path(basedir + '/OUTCAR')
    qe = Path(basedir + '/run_aseQE.py')

    if vasp.is_file():
        # VASP calculation
        ## Figure out if it is a implicit / explicit calc
        homedir = basedir #+ '/0' + str(images) + '/' 
        test = subprocess.Popen("grep 'LSOL' "+homedir+"/OUTCAR",\
                shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        if test == '':
            implicit = False
            # Getting WF for explicit case
            try:
                wf = vasp_wf.get_wf_explicit(homedir)
            except FileNotFoundError:
                wf = np.nan
            except IndexError:
                wf = np.nan
        elif test.split()[2] == 'T':
            implicit = True
            # Getting WF for implicit case
            wf = vasp_wf.get_wf_implicit(homedir, neb=True)
        print(wf)
    elif qe.is_file():
        homedir = basedir + 'out_000' + str(image)
        # Quantum espresso calculation
        test = subprocess.Popen("grep 'Environ Module' "+homedir+"/log",\
                shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        if test == '':
            implicit = False
            try:
                wf = qe_wf.get_wf_explicit(homedir)
            except FileNotFoundError:
                wf = np.nan
        else:
            implicit = True
            wf = qe_wf.get_wf_implicit(homedir)

    return [implicit, wf]


def get_charge(homedir):
    qe = Path(homedir + '/run_aseQE.py')
    vasp = Path(homedir + '/OUTCAR')

    if vasp.is_file():
        test = subprocess.Popen("grep 'NELECT' "+homedir+"/OUTCAR",\
                shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        if test == '':
            charge = 0
        else:
            split = test.split('=')[1].split(',')[0].split()[0]
            charge_text = split
            charge = float(charge_text.replace('d', 'e'))
    elif qe.is_file():
        # Returns tot_charge based on pw input
        # If it doesnt find it return 0 ie uncharged
        outdir = homedir + 'out_0000/' 
        # Looks for tot_charge flag in pw input
        f = open(outdir + 'pw.inp')
        lines = f.readlines()
        for line in lines:
            if 'tot_charge' in line:
                split = (line.split('='))[1].split(',')
                charge_text = split[0]
                charge = float(charge_text.replace('d', 'e'))
        f.close()
    if charge:
        return charge
    else:
        return 0.0


def get_other_info(homedir):
    
    # Getting pw cutoff
    # functional used
    # Pseudopotential / PAW potential used
    
    qe = Path(homedir + '/run_aseQE.py')
    vasp = Path(homedir + 'OUTCAR')

    if vasp.is_file():
        test = subprocess.Popen("grep 'ENCUT' "+homedir+"OUTCAR",\
                shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        pw = float(test.split()[2])
        test = subprocess.Popen("grep 'GGA' "+homedir+"OUTCAR",\
                shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        functional = test.split()[-1]
        if functional == '--' or functional == 'PE':
            lhfcalc = subprocess.Popen("grep 'LHFCALC' "+homedir+"OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            hfscreen = subprocess.Popen("grep 'HFSCREEN=' "+homedir+"OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if lhfcalc.split()[2] == 'T' and hfscreen.split()[1] == '0.2000':
                functional = 'HSE'
        test = subprocess.Popen("grep ' LDAUTYPE' "+homedir+"OUTCAR",\
                shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        if test == '':
            ldautype = np.nan
            angular_moment = np.nan
            ldauU = 0
            ldauJ = np.nan
            ldau_data = {'ldautype':ldautype, 'angular_moment':angular_moment, \
                    'U':ldauU, 'J':ldauJ}
        else:
            # GGA U is active
            ldautype = int(test.split()[-1])
            test = subprocess.Popen("grep 'angular momentum for each species LDAUL' "+homedir+"/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            temp_string = test.split('=')[-1].split()
            angular_moment = [ float(i) for i in temp_string ]
            test = subprocess.Popen("grep 'LDAUU' "+homedir+"/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            temp_string = test.split('=')[-1].split()
            ldauU = [ float(i) for i in temp_string ]
            test = subprocess.Popen("grep 'LDAUJ' "+homedir+"/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            temp_string = test.split('=')[-1].split()
            ldauJ = [ float(i) for i in temp_string ]

            ldau_data = {'ldautype':ldautype, 'angular_moment':angular_moment, \
                    'U':ldauU, 'J':ldauJ}

        paw = 'vasp'

        data = {'pw_cutoff':pw, 'functional':functional, 'paw':paw,\
                'ldau_data':ldau_data}

    elif qe.is_file():
        data = {}
        paw = []
        functional = ''
        outdir = homedir + 'out_0000/' 

        f = open(outdir + 'pw.inp')
        lines = f.readlines()
        for line in lines:
            if 'ecutwfc' in line:
                split = (line.split('='))[1].split(',')
                text = split[0]
                pw_convert = float(text.replace('d', 'e')) * units.Ry
                pw = int(pw_convert)
            if 'input_dft' in line:
                split = (line.split('='))[1].split(',')
                functional = str(split[0].split("'")[1])
            if 'pseudo_dir' in line:
                split = (line.split('='))[1].split(',')[0].split('/')
                paw = split[-1].split("'")[0]
        data = {'pw_cutoff':pw, 'functional':functional, 'paw':paw}

    return data 


def get_dirs(homedir):
    # Given a path, returns a list of directories
    # Also returns if state information is provided

    send_direct = []

    for direct in os.listdir(homedir):
        if os.path.isdir(homedir + '/' + direct) and direct != 'archive' and \
                direct != 'intermediates' and direct != 'starting_structures':
            send_direct.append(homedir + '/' +  str(direct) )
    if os.path.exists(homedir + 'COa_COg') or os.path.exists(homedir + 'CO2g_COa'):
        present_states = True
    else:
        present_states = False
    return [send_direct, present_states]
    
    
def check_climb(homedir):
    qe = Path(homedir + '/run_aseQE.py')
    vasp = Path(homedir + '..' +  '/INCAR')

    if vasp.is_file():
        test = subprocess.Popen("grep 'LCLIMB' "+homedir+"../INCAR",\
                shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
        if test.split() == '':
            climb = False
        elif test.split()[-1] == '.TRUE.':
            climb=True
        elif test.split()[-1] == '.FALSE.':
            climb=False

    elif qe.is_file():
       # Check if neb with clim turned on or not
       # read in all slurm files
       fslurm = glob(homedir + '/slurm*')
       
       f = open(Path(fslurm[-1]))

       lines = f.readlines()
       for line in lines:
           line = line.rstrip()
           if "climb mode on" in line:
               climb = True
           elif "no climbing" in line:
               climb = False
    return climb

def number_of_images(homedir):
    qe = Path(homedir + '/run_aseQE.py')
    vasp = Path(homedir + 'INCAR')
    # Find the number of images based on the number of neb files
    if vasp.is_file():
        nebfiles = glob(homedir +  '0*/')
        images = len(nebfiles)
    elif qe.is_file():
        nebfiles = glob(homedir + 'neb*.traj')
        images = len(nebfiles)
    return images

        
if __name__ == '__main__':
    """
    go through directories and get atoms objects for each neb
    parse through the atoms object and get all relevant quantitites
    """
    dbname = str(sys.argv[1])
    levels_max = int(sys.argv[2])

    levels = np.arange(levels_max)

    # Getting directories to parse
    homedirs = []

    # All path information
    paths = {}
    # Does the path contain states info
    flag_states = []
    """ 
    get all homedirs 
    TODO: Change this at some point
    """
    # Goes through each file level and collects path information
    for i in range(len(levels)):
        level = levels[i] 
        print(paths, level)
        if level == 0:
            current_list, flag_s = get_dirs('.')
            paths[level] = current_list
            flag_states.append(flag_s)
        elif level > 0:
            accumulated_list = []
            for j in range(len(paths[level-1])):
                current_list, flag_s = get_dirs(paths[level-1][j])
                accumulated_list.append(current_list)
            flag_states.append(flag_s)
            saved_list = [item for sublist in accumulated_list for item in sublist]
            paths[level] = saved_list
    # Feeding all paths into homedir
    homedirs = paths[levels_max-1] 
    # Get states
    states = []
    db = connect(dbname)

    # Loop through all final paths and add into to the database
    for i in range(len(homedirs)):
        homedir = homedirs[i] + '/' 
        print("directories parsed:")
        print(homedir)
        # For a give run, find things common for all images
        charge = get_charge(homedir)
        climb = check_climb(homedir)
        other_info = get_other_info(homedir)
        functional = other_info['functional']
        pw = other_info['pw_cutoff']
        paw = other_info['paw']
        # Now go for things that are specific for each image
        ## find the number of images
        image = homedir.split('/')[-2]#number_of_images(homedir + '/../')
        print('image number: %s'%image)
        #    image = i 
        atoms_object = get_atoms(homedir, image)
        if atoms_object == "":
            print("Something wrong here")
        else:
            implicit, wf = get_wf_parsed(homedir, image)
        # Postion of the states variable
        split_list = homedir.split('/')
        for lst in split_list:
            if 'slab' in lst or 'CO' in lst or 'CH' in lst or 'H' in lst or 'Li' in lst:
                states_names = lst
            if 'DV' in lst or 'SV' in lst or 'SW' in lst or '11' in lst or 'Co' in lst:
                facet_names = lst
            if 'x' in lst: # change this when you can
                cell_size = lst
        states.append(states_names)
        db.write(atoms_object, \
                implicit=implicit,  wf=wf,\
                tot_charge=float(charge),\
                states=states_names, \
                functional = functional,\
                pw_cutoff = float(pw), \
                paw=paw,\
                image=float(image), climb=climb,
                data={'other_info':other_info,\
            }
        )


