#!/usr/bin/python


""" Functions for helping restart calculations """
import sys, subprocess, os
import numpy as np
from pprint import pprint
from ase.io import read
from pathlib import Path

def get_vasp_nelect0(atoms):
    import pickle
    import numpy as np
    # Given an atoms object it will return the nelect at zero charge
    # Based on default vasp potentials
    with open('/Users/vijays/Documents/tools/scripts/nelect0.pickle','rb') as handle:
        nelect0 = pickle.load(handle)
    default_nelect = []
    for i in range(len(atoms)):
        default_nelect.append(nelect0[atoms[i].symbol])
    default_nelect = np.array(default_nelect)
    return default_nelect.sum()

def _check_code(homedir):
    # Checks which code was used - there are only two options
    codes = {
    'qe':[ 'qn.traj', 'bfgs.traj', 'spe.traj', 'log'],
    'vasp':['OUTCAR', 'vasprun.xml', 'CONTCAR', 'POSCAR'],
    #'gpaw':['qn_gpaw.traj', 'init_gpaw.traj'],
        }
    for code in codes:
        for f in codes[code]:
            filetype = Path(homedir + '/' + f)
            if filetype.is_file():
                print('This is a ' + code + ' calculation')
                code_used = code
                break

    return code_used

def _status_calc(homedir, user, pwd):
    # Check if calculation is queued or running

    current_jobs = subprocess.check_output(['squeue','-u','%s'%user,'-o','%.18i %.9P %.8j %.8u %.2t %.10M %.6D %.20V %.20S %Z']).decode('utf-    8').split('\n')[:-1]
    current_jobs = [job for job in current_jobs if job != None]
    id_all = current_jobs[:18]
    if '_' in id_all:
        id_all = id_all.split('_')[0]
    current_ids = []
    for current_job in current_jobs[1:]:
        ids = current_job.split()[0]
        if '_' not in ids: # cant handle slurm arrays
            try:
                current_ids.append(str(int(ids)))
            except TypeError:
                continue
    current_dirs = []
    for job in current_jobs:
        current_dirs.append(job.split()[9].strip())
    full_homedir = os.path.join(pwd, homedir)
    if os.path.normpath(full_homedir) in current_dirs:
        status = 'queue'
    else:
        # Check which code is used 
        code = _check_code(homedir)
        if code == 'vasp':
            # Check if the calculation has completed
            test = subprocess.Popen(\
                    "grep 'reached required accuracy - stopping structural energy minimisation' "\
                    +homedir+"/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if bool(test):
                status = 'completed'
            else:
                # Check if it was a single point
                test = subprocess.Popen(\
                        "grep 'NSW' "\
                        +homedir+"/INCAR",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                try:
                    if int(test.split(' ')[-1]) > 0:
                        status = 'failed'
                    else:
                        status = 'completed'
                except ValueError:
                    status = 'failed'

        elif code == 'qe':
            # If it is an optimization check if forces are below 0.025 eV / A
            if Path(homedir + '/qn.traj').exists():
                atoms = read(homedir + '/qn.traj')
                if max([np.linalg.norm(a) for a in atoms.get_forces()]) < 0.025:
                    status = 'completed'
                else:
                    status = 'failed'
            else:
                # Check if it is a single point 
                if Path(homedir + '/spe.traj').exists():
                    status = 'completed'
                else:
                    status = 'failed'

    return status

def _restart_calc(homedir):

    code = _check_code(homedir)

    if code == 'vasp':
        # Restart by creating an archive directory
        os.system('mkdir -p ' + homedir + '/archive')
        try:
            atoms = read(homedir + '/CONTCAR')
        except:
            atoms = read(homedir + '/POSCAR')
        os.system('mv ' + homedir + '/** ' + homedir + '/archive/')

        # Keep the WAVECAR run_vasp.py and submit script in the directory
        os.system('mv ' + homedir + '/archive/WAVECAR ' + homedir)
        os.system('cp ' + homedir + '/archive/{run_vasp.py,submit_scriptVASP.sh} ' + homedir)

        atoms.write(homedir + '/init.traj')
        atoms.write(homedir + '/CONTCAR') # sometimes scripts use CONTCAR

        # write it out to ensure easy restarts 
        with open('run_list_fail', 'a') as f:
            f.write(homedir)
            f.write('\n')

        owd = os.getcwd()
        os.chdir(homedir)
        os.system('zip -r archive.zip archive/')
        os.system('sbatch submit_scriptVASP.sh')
        #os.system('rm -r archive')
        os.chdir(owd)

    elif code == 'qe':
        # Restart by creating an archive directory
        os.system('mkdir -p ' + homedir + '/archive')
        try:
            atoms = read(homedir + '/qn.traj')
        except:
            atoms = read(homedir + '/init.traj')
        os.system('mv ' + homedir + '/** ' + homedir + '/archive/')

        os.system('cp ' + homedir + '/archive/{run_aseQE.py,submit_script*} ' + homedir)

        atoms.write(homedir + '/init.traj')

        with open('run_list_fail', 'a') as f:
            f.write(homedir)
            f.write('\n')

        owd = os.getcwd()
        os.chdir(homedir)
        os.system('zip -r archive.zip archive/')
        os.system('sbatch submit_scriptQE61.sh')
        os.system('rm -r archive')
        os.chdir(owd)

