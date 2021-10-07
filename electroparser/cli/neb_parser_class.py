#!/usr/bin/python


""" Parser class that takes care of all DFT calculations """


from ase.io import read, write
from ase.db import connect
import numpy as np
import os, sys
from pathlib import Path
from electroparser.cli.parse_wf import get_wf as qe_wf
from electroparser.cli.parse_wf import get_wf_initial_implicit
from electroparser.cli import parse_wf_vasp as vasp_wf
import re
from ase import units
import glob
import subprocess
from ase.vibrations import Vibrations
from glob import glob
from electroparser.cli.useful_functions import d_band_info, get_bands_DOS

class NebParser:
    # Class that parses a neb
    def __init__(self, homedir):
        ### INPUTS
        self.homedir = homedir # the home directory to parse form

        ### OUTPUTS
        self.directories = []
        self.atoms = [] # atoms object that starts empty append images as we go
        self.code = '' # code used
        self.implicit = False # Check if implicit calculation was performed
        self.wf = [] # Progression of work function of calculation
        self.field = 0.0 # Electric field in calculation (same for neb)
        self.dipole = [] # Dipole moment from DFT across images
        self.charge = 0.0 # Excess charge added to calc
        self.input_file = [] # Contents of input file
        self.paw = 0 # Plane wave cutoff
        self.functional = '' # Functional used
        self.no_images = 0 # number of images
        self.ldau_data = {} # All U related data

        ### UNITS
        self.units = {'debye2eA':0.20819434}

        ### Execute all functions
        self.get_atoms()
        self.get_wf()
        self.get_field()
        self.get_dipole()
        self.get_charge()
        self.get_input_file()
        self.get_other_info()


    # Get the atoms object from whatever files are in the homedir
    def get_atoms(self):

        # Written in order of which it should parse
        codes = {
        'qe':[ 'neb.traj', 'neb0.traj'],
        'vasp':['00'],
        #'gpaw':['qn_gpaw.traj', 'init_gpaw.traj'],
            }

        for code in codes:
            for f in codes[code]:
                filetype = Path(self.homedir + '/' +  f)
                if filetype.is_file() or filetype.is_dir():
                    print('Storing ' + f + ' as atoms object')
                    self.code = code
                    break

        if self.code == 'vasp':
            image_dirs = glob(self.homedir + '0?/')
            image_dirs = sorted(image_dirs)
            try:
                for image in image_dirs:
                    self.atoms.append(read(image +  '/OUTCAR'))
            except ValueError:
                print('WARNING: could not store ' + f )


        elif self.code == 'qe':
            image_dirs = glob(self.homedir + 'neb?.traj')
            print(image_dirs)
            image_dirs = sorted(image_dirs)
            for image in image_dirs:
                self.atoms.append(read(image))

        # Number of images for the purpose of iterating
        self.no_images = len(image_dirs)


    # Get the workfunction if files are available
    def get_wf(self):
        if self.code == 'vasp':
            print(self.homedir)
            test = subprocess.Popen("grep 'LSOL' "+self.homedir+"01/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.implicit = False

                # Getting WF for explicit case
                for index in range(self.no_images):
                    if index < 10:
                        image_no = '0' + str(index)
                    else:
                        image_no = str(image_no)
                    try:
                        self.wf.append(vasp_wf.get_wf_explicit(self.homedir+image_no+'/'))
                    except (IndexError, FileNotFoundError):
                        self.wf.append(np.nan)
            elif test.split()[2] == 'T':
                self.implicit = True
                # Getting WF for implicit case
                self.wf = vasp_wf.get_wf_implicit(self.homedir)

        elif self.code == 'qe':
            test = subprocess.Popen("grep 'Environ Module' "+self.homedir+"/log",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.implicit = False
                for index in range(self.no_images):
                    try:
                        self.wf.append(qe_wf.get_wf_explicit(self.homedir))
                    except FileNotFoundError:
                        self.wf.append(np.nan)
            else:
                self.implicit = True
                try:
                    for index in range(self.no_images):
                        self.wf.append(qe_wf.get_wf_implicit(self.homedir))
                except FileNotFoundError:
                    self.wf.append(np.nan)



        elif self.code == 'gpaw':
            test = subprocess.Popen("grep 'Dipole-l' "+self.homedir+"/gpawlog",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.implicit = False
                self.wf = np.nan
            else:
                self.implicit = True
                self.wf = float(test.split()[-2])

    # Gets the electric field in the calculation
    def get_field(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'EFIELD' "+self.homedir+"/01/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                # No field present
                self.field = 0.
            else:
                # There is a field
                self.field = float(test.split()[2])

        elif self.code == 'qe':
                test = subprocess.Popen("grep 'eamp' "+self.homedir+"/pw.inp",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                if test == '':
                    self.field = 0
                else:
                    split = test.split('=')[1].split(',')
                    field_text = split[0]
                    self.field = float(field_text.replace('d', 'e'))

    # Gets the dipole moment in the calculation
    def get_dipole(self):
        if self.code == 'vasp':
            print(self.no_images)
            for index in range(self.no_images):
                if index < 10:
                    image_no = '0' + str(index)
                else:
                    image_no = str(image_no)

                test = subprocess.Popen("grep 'dipolmoment' "+self.homedir+'/'+image_no+"/OUTCAR | tail -1" ,\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')

                if test == '':
                    # No dipole from the field
                    self.dipole.append(np.nan)
                else:
                    self.dipole.append(float(test.split()[-4]))

        elif self.code == 'qe':
            test = subprocess.Popen("grep -A1 'dipole' "+self.homedir+"/log",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                # No dipole from the field
                self.dipole = np.nan * np.ones(self.no_images)
            else:
                # TODO: Find a better way of parsing the dipole moment
                self.dipole = float(test.split()[-2]) * self.units['debye2eA']

    def get_charge(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'NELECT' "+self.homedir+"/01/OUTCAR | tail -1" ,\
            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                # No dipole from the field
                self.charge = 0.0
            else:
                self.charge = float(test.split()[2])
        if self.code == 'qe':
            test = subprocess.Popen("grep 'tot_charge' "+self.homedir+"/pw.inp",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.charge = 0
            else:
                split = test.split('=')[1].split(',')
                charge_text = split[0]
                self.charge = float(charge_text.replace('d', 'e'))

    def get_input_file(self):
        if self.code == 'vasp':
            try:
                f = open(self.homedir + 'INCAR')
                self.input_file = f.readlines()
            except:
                print('WARNING: no input file stored')
        elif self.code == 'qe':
            try:
                f1 = open(self.homedir + 'pw.inp')
            except FileNotFoundError:
                self.input_file = []

    def get_other_info(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'ENCUT' "+self.homedir+"/01/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test != '':
                self.pw = float(test.split()[2])
                test = subprocess.Popen("grep 'GGA type' "+self.homedir+"/01/OUTCAR",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                self.functional = test.split()[2]
                # check if vdw is active
                test = subprocess.Popen("grep 'IVDW' "+self.homedir+"/01/OUTCAR",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                if test != '':
                    dispersion_type = test.split()[2]
                    if dispersion_type == '11':
                        self.functional = functional + '+' + 'D3'
                    elif dispersion_type == '12':
                        self.functional = functional + '+' + 'D3BJ'

                if self.functional == '--' or self.functional == 'PE':
                    lhfcalc = subprocess.Popen("grep 'LHFCALC' "+self.homedir+"/01/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    hfscreen = subprocess.Popen("grep 'HFSCREEN=' "+self.homedir+"/01/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    if lhfcalc.split()[2] == 'T' and hfscreen.split()[1] == '0.2000':
                        self.functional = 'HSE'
                test = subprocess.Popen("grep ' LDAUTYPE' "+self.homedir+"/01/OUTCAR",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                if test == '':
                    ldautype = np.nan
                    angular_moment = np.nan
                    ldauU = 0
                    ldauJ = np.nan
                    self.ldau_data = {'ldautype':ldautype, 'angular_moment':angular_moment, \
                            'U':ldauU, 'J':ldauJ}
                else:
                    # GGA U is active
                    ldautype = int(test.split()[-1])
                    test = subprocess.Popen("grep 'angular momentum for each species LDAUL' "+self.homedir+"/01/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    temp_string = test.split('=')[-1].split()
                    angular_moment = [ float(i) for i in temp_string ]
                    test = subprocess.Popen("grep 'LDAUU' "+self.homedir+"/01/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    temp_string = test.split('=')[-1].split()
                    ldauU = [ float(i) for i in temp_string ]
                    test = subprocess.Popen("grep 'LDAUJ' "+self.homedir+"/01/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    temp_string = test.split('=')[-1].split()
                    ldauJ = [ float(i) for i in temp_string ]

                    self.ldau_data = {'ldautype':ldautype, 'angular_moment':angular_moment, \
                            'U':ldauU, 'J':ldauJ}

            self.paw = 'vasp'

        elif self.code == 'qe':
            try:
                f = open(self.homedir + 'pw.inp')
                lines = f.readlines()
                for line in lines:
                    if 'ecutwfc' in line:
                        split = (line.split('='))[1].split(',')
                        text = split[0]
                        pw_convert = float(text.replace('d', 'e')) * units.Ry
                        self.pw = int(pw_convert)
                    if 'input_dft' in line:
                        split = (line.split('='))[1].split(',')
                        self.functional = str(split[0].split("'")[1])
                    if 'pseudo_dir' in line:
                        split = (line.split('='))[1].split(',')[0].split('/')
                        self.paw = split[-1].split("'")[0]
            except FileNotFoundError:
                self.paw = 'qe_unknown_paw'
                self.functional = 'no_input_found'
                self.pw = np.nan
