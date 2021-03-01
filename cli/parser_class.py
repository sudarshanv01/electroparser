#!/usr/bin/python


""" Parser class that takes care of all DFT calculations """


from ase.io import read, write
from ase.db import connect
import numpy as np
import os, sys
from pathlib import Path
sys.path.append('/home/cat/vijays/tools/scripts')
import get_wf as qe_wf
from get_wf import get_wf_initial_implicit
import get_wf_vasp as vasp_wf
import re
from ase import units
import glob
from ase.io.bader import attach_charges
import subprocess
from ase.vibrations import Vibrations
from glob import glob
from useful_functions import d_band_info, get_bands_DOS
from copy import deepcopy 
import json
import ase

class Parser:

    def __init__(self, homedir):
        ### INPUTS
        self.homedir = homedir # the home directory to parse form

        ### OUTPUTS
        self.directories = []
        self.atoms = [] # atoms object that starts empty
        self.code = '' # code used
        self.implicit = False # Check if implicit calculation was performed
        self.wf = 0.0 # Work function of calculation
        self.field = 0.0 # Electric field in calculation
        self.dipole = 0.0 # Dipole moment
        self.charge = 0.0 # Excess charge added to calc
        self.input_file = [] # Contents of input file
        self.paw = 0 # Plane wave cutoff
        self.functional = '' # Functional used
        self.ldau_data = {} # All U related data
        self.energy = [] # Only in the case that POSCAR is stored
        self.spin_pol = True # Is this a spin polarized calculation
        self.ldau = False # Is this a LDA+U calculation
        self.d_centre = 0.0 # d band centre
        self.pw = 0.0 # plane wae cutoff
        self.d_centre_up = 0.0 # d band centre up
        self.d_centre_down = 0.0 # d band centre down
        self.width = 0.0 # width of d band
        self.width_up = 0.0 # width of d band
        self.width_down = 0.0 # width of d band
        self.max_hilbert = 0.0 #max of hilbert
        self.max_hilbert_up = 0.0 # max hilbert up
        self.max_hilbert_down = 0.0 # max hilbert down
        self.occupancy_up = 0.0 # Occupancy of spin up state
        self.occupancy_down = 0.0 # Occupancy of spin down state
        self.occupancy = 0.0 # Total occupant for non spin polarized cases
        self.vibrations = [] # List of vibrational frequencies if any
        self.bader_atoms = [] # atoms object into which bader charges go
        self.bader = False # Bool to indicate that the bader charges exits
        self.ir_intensity = 0.0
        self.raman_intensity = 0.0
        self.eigenval = {} #eigenvalues if available
        self.pdos = {} # projected dos if available
        self.Vxy = {} # xy-averages potential if available
        self.sp_energy = 0.0 # Energy of sp if relaxation

        ### UNITS
        self.units = {'debye2eA':0.20819434}

        ### Execute all functions
        self.get_atoms()
        self.get_wf()
        self.get_d_info()
        self.get_field()
        self.get_dipole()
        self.get_charge()
        self.get_input_file()
        self.get_other_info()
        self.get_vibrations()
        #self.get_spectra() # Gets the IR and Raman spectra if exists
        self.attach_bader()
        self.get_eigenvals()
        self.get_pdos()
        self.get_Vxy()



    # Get the atoms object from whatever files are in the homedir
    def get_atoms(self):

        # Written in order of which it should parse
        codes = {
        'qe':[ 'qn.traj', 'bfgs.traj', 'spe.traj', 'aiida.out', 'log'],
        'vasp':['OUTCAR', 'vasprun.xml', 'CONTCAR', 'POSCAR'],
        'gpaw':['gpaw.traj', 'out.txt'],
        #'gpaw':['qn_gpaw.traj', 'init_gpaw.traj'],
            }
        success = False
        for code in codes:
            if success:
                break
            for f in codes[code]:
                filetype = Path(self.homedir + f)
                if filetype.is_file():
                    success = True
                    print('Storing ' + f + ' as atoms object')
                    try:
                        self.atoms = read(filetype)
                    except (ValueError, StopIteration):
                        print('WARNING: could not store ' + f )
                    self.code = code
                    if f == 'POSCAR':
                        # Need to store the energy here
                        with open(self.homedir + 'energy.out') as ener_f:
                            s = ener_f.read()
                            self.energy = float(s)
                    break
        ## check if information about the first run is available 
        try:
            all_atoms = read(filetype, ':')
            self.sp_energy = all_atoms[0].get_potential_energy()
        except ase.calculators.calculator.PropertyNotImplementedError:
            pass
        except FileNotFoundError:
            pass
    
    # Get the d band center and d band width
    def get_d_info(self):
        try:
            with open(self.homedir + 'd_center.txt') as f:
                d_centres = f.readlines()
                if len(d_centres) == 2:
                    self.spin_pol = True
                    self.d_centre_up = float(d_centres[0])
                    self.d_centre_down = float(d_centres[1])
                else:
                    self.spin_pol = False
                    self.d_centre = float(d_centres[0])
        except:
            self.spin_pol = np.nan
        try:
            with open(self.homedir + 'max_hilbert.txt') as f:
                max_hils = f.readlines()
                if len(max_hils) == 2:
                    self.spin_pol = True
                    self.max_hilbert_up = float(max_hils[0])
                    self.max_hilbert_down = float(max_hils[1])
                else:
                    self.max_hilbert = float(max_hils[0].split('\n')[0])
        except FileNotFoundError:
            # print('No information about hilbert maximum')
            self.spin_pol = np.nan

        try:
            with open(self.homedir + 'occupancy.txt') as f:
                occupancies = f.readlines()
                if len(occupancies) == 2:
                    self.occupancy_up = abs(float(occupancies[0]))
                    self.occupancy_down = abs(float(occupancies[1]))
                else:
                    self.occupancy = abs(float(occupancies[0].split('\n')[0]))
        except FileNotFoundError:
            # print('No occupancy data for spinpol calculation')
            pass
            # self.spin_pol = np.nan
        
        try: 
            with open(self.homedir + 'd_width.txt') as f:
                width = f.readlines()
                if len(width) == 2:
                    self.width_up = abs(float(width[0]))
                    self.width_down = abs(float(width[1]))
                else:
                    self.width = abs(float(width[0].split('\n')[0]))
        except FileNotFoundError:
            # print('No width data')
            pass
            # self.spin_pol = np.nan

    def get_eigenvals(self):
        ## get the eigenvalues for a GPAW calculation
        if self.code == 'gpaw':
            files = glob(os.path.join(self.homedir, 'eigenvals*'))
            data = {}
            for f in files:
                filename = f.split('/')[-1]
                data[filename] = np.loadtxt(f)
            self.eigenval = data

    def get_pdos(self):
        if 'pdos.json' in os.listdir(self.homedir):
            with open(os.path.join(self.homedir, 'pdos.json'), 'r') as handle:
                self.pdos = json.load(handle)

    def get_Vxy(self):
        if 'Vxy.json' in os.listdir(self.homedir):
            with open(os.path.join(self.homedir, 'Vxy.json'), 'r') as handle:
                self.Vxy = json.load(handle)

    # Get the workfunction if files are available
    def get_wf(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'LSOL' "+self.homedir+"/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.implicit = False
                # Getting WF for explicit case
                try:
                    self.wf = vasp_wf.get_wf_explicit(self.homedir)
                except FileNotFoundError:
                    self.wf = np.nan
                except IndexError:
                    self.wf = np.nan
                except ValueError:
                    self.wf = np.nan
            elif test.split()[2] == 'T':
                self.implicit = True
                # Getting WF for implicit case
                self.wf = vasp_wf.get_wf_implicit(self.homedir)

        elif self.code == 'qe':
            test = subprocess.Popen("grep 'Environ Module' "+self.homedir+"/log",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.implicit = False
                try:
                    self.wf = qe_wf.get_wf_explicit(self.homedir)
                except (FileNotFoundError,UnboundLocalError):
                    self.wf = np.nan
            else:
                self.implicit = True
                self.wf = qe_wf.get_wf_implicit(self.homedir)

        elif self.code == 'gpaw':
            possible_gpaw_filename = ['gpaw.txt', 'gpaw.log', 'out.txt']
            for f in possible_gpaw_filename:
                try:
                    subprocess.check_output(['grep',  "Dipole-layer",  self.homedir+f])
                    filename = f
                    break
                except subprocess.CalledProcessError:
                    continue
            try:
                test = subprocess.Popen("grep 'Dipole-layer corrected work functions' "+self.homedir+"/"+filename,\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            except UnboundLocalError:
                test = ''
            if test == '':
                self.wf = np.nan
            else:
                self.wf = float(test.split()[-2])
            try:
                test = subprocess.Popen("grep 'SJM' "+self.homedir+"/"+filename,\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            except UnboundLocalError:
                test = ''
            if test == '':
                self.implicit = False
            else:
                self.implicit = True

    # Gets the electric field in the calculation
    def get_field(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'EFIELD' "+self.homedir+"/OUTCAR",\
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
                
        elif self.code == 'gpaw':
            possible_gpaw_filename = ['gpaw.txt', 'gpaw.log', 'out.txt']
            for f in possible_gpaw_filename:
                try:
                    subprocess.check_output(['grep',  "Dipole-layer",  self.homedir+f])
                    filename = f
                    break
                except subprocess.CalledProcessError:
                    continue
            try:
                test = subprocess.Popen("grep ' Constant electric field:' "+self.homedir+"/"+filename,\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            except UnboundLocalError:
                test = ''
            if test == '':
                self.field = 0.0
            else:
                self.field = float(test.split(')')[0].split()[-1])

    # Gets the dipole moment in the calculation
    def get_dipole(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'dipolmoment' "+self.homedir+"/OUTCAR | tail -1" ,\
            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                # No dipole from the field
                self.dipole = np.nan
            else:
                self.dipole = float(test.split()[-4])

        elif self.code == 'qe':
            test = subprocess.Popen("grep -A1 'dipole' "+self.homedir+"/log",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                # No dipole from the field
                self.dipole = np.nan
            else:
                # TODO: Find a better way of parsing the dipole moment
                self.dipole = float(test.split()[-2]) * self.units['debye2eA']
        elif self.code == 'gpaw':
            possible_gpaw_filename = ['gpaw.txt', 'gpaw.log', 'out.txt']
            for f in possible_gpaw_filename:
                try:
                    subprocess.check_output(['grep',  "Dipole-layer",  self.homedir+f])
                    filename = f
                    break
                except subprocess.CalledProcessError:
                    continue
            try:
                test = subprocess.Popen("grep 'Dipole moment' "+self.homedir+"/"+filename,\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                self.dipole = float(test.split()[-2].replace(')',''))
            except UnboundLocalError:
                pass

    def get_charge(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'NELECT' "+self.homedir+"/OUTCAR | tail -1" ,\
            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                # No dipole from the field
                self.charge = 0.0
            else:
                self.charge = float(test.split()[2])
        elif self.code == 'qe':
            test = subprocess.Popen("grep 'tot_charge' "+self.homedir+"/pw.inp",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.charge = 0
            else:
                split = test.split('=')[1].split(',')
                charge_text = split[0]
                self.charge = float(charge_text.replace('d', 'e'))
                
        elif self.code == 'gpaw':
            try:
                subprocess.check_output(['grep',  "Dipole-layer",  self.homedir+"gpaw.txt"])
                filename = 'gpaw.txt'
            except subprocess.CalledProcessError:
                try:
                    subprocess.check_output(['grep',  "Dipole-layer",  self.homedir+'gpaw.log'])
                    filename = 'gpaw.log'
                except subprocess.CalledProcessError:
                    #filename = 'gpaw.log'
                    filename = 'log'

            test = subprocess.Popen("grep 'Current number of Excess Electrons' "+self.homedir+"/"+filename,\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test == '':
                self.charge = 0
            else:
                self.charge = float(test.split()[-1])

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
                f1 = open(self.homedir + 'aiida.in')
            try:
                f2 = open(self.homedir + 'environ.in')
                self.input_file = f1.readlines() + f2.readlines()
            except FileNotFoundError:
                self.input_file = f1.readlines()

    def get_other_info(self):
        if self.code == 'vasp':
            test = subprocess.Popen("grep 'ENCUT' "+self.homedir+"/OUTCAR",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test != '':
                self.pw = float(test.split()[2])
                test = subprocess.Popen("grep 'GGA type' "+self.homedir+"/OUTCAR",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                self.functional = test.split()[2]
                # check if vdw is active
                test = subprocess.Popen("grep 'IVDW' "+self.homedir+"/OUTCAR",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                if test != '':
                    dispersion_type = test.split()[2]
                    if dispersion_type == '11':
                        self.functional = self.functional + '+' + 'D3'
                    elif dispersion_type == '12':
                        self.functional = self.functional + '+' + 'D3BJ'

                if self.functional == '--' or self.functional == 'PE':
                    lhfcalc = subprocess.Popen("grep 'LHFCALC' "+self.homedir+"/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    hfscreen = subprocess.Popen("grep 'HFSCREEN=' "+self.homedir+"/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    if lhfcalc.split()[2] == 'T' and hfscreen.split()[1] == '0.2000':
                        self.functional = 'HSE'
                test = subprocess.Popen("grep ' LDAUTYPE' "+self.homedir+"/OUTCAR",\
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
                    self.ldau = True
                    ldautype = int(test.split()[-1])
                    test = subprocess.Popen("grep 'angular momentum for each species LDAUL' "+self.homedir+"/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    temp_string = test.split('=')[-1].split()
                    angular_moment = [ float(i) for i in temp_string ]
                    test = subprocess.Popen("grep 'LDAUU' "+self.homedir+"/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    temp_string = test.split('=')[-1].split()
                    ldauU = [ float(i) for i in temp_string ]
                    test = subprocess.Popen("grep 'LDAUJ' "+self.homedir+"/OUTCAR",\
                            shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                    temp_string = test.split('=')[-1].split()
                    ldauJ = [ float(i) for i in temp_string ]

                    self.ldau_data = {'ldautype':ldautype, 'angular_moment':angular_moment, \
                            'U':ldauU, 'J':ldauJ}

            self.paw = 'vasp'

        elif self.code == 'qe':
            try:
                f = open(self.homedir + 'pw.inp')
            except FileNotFoundError:
                f = open(self.homedir + 'aiida.in')
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

            # Check if there is a Hubbard U correction used 
            test = subprocess.Popen("grep 'lda_plus_u' "+self.homedir+"/pw.inp",\
                    shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
            if test != '':
                # A Hubbard U calculation was run
                self.ldau = True
                test = subprocess.Popen("grep 'Hubbard_U' "+self.homedir+"/pw.inp",\
                        shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                atom_U_list = test.split('\n')
                ldauU = [ float(a.split('=')[-1].replace('d', 'e').replace(',','')) for a in test.split('\n')[0:-1] ]
                self.ldau_data = {'U':ldauU}

    def get_vibrations(self):

        if self.code == 'vasp':
            if Path(self.homedir + '/vasp_freq').is_dir():
                # Needs a vib directory where VASP vibration calc
                # has run 
                test = subprocess.Popen("grep 'cm-1' "+self.homedir+"/vasp_freq/OUTCAR",\
                                 shell=True,stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
                ltest = test.split('\n')
                freq = []
                im_freq = []
                for line in ltest:
                    data = line.split()
                    if 'f/i=' not in data and len(data) > 0:
                        # real frequencies
                        freq.append(float(data[-4]))
                    elif 'f/i=' in data and len(data) > 0:
                        im_freq.append(float(data[-4]))
                # Store only real vibrations
                self.vibrations = freq
    
    def get_spectra(self):
        # Gets the raman and IR intensities 
        with open(self.homedir + '/ir_intensity.txt') as f:
            self.ir_intensity = float(f.readline())
        with open(self.homedir + '/raman_intensity.txt') as f:
            self.raman_intensity = float(f.readline())

    def attach_bader(self):
        acf= Path(self.homedir + '/ACF.dat')
        acf_folder = Path(self.homedir + '/bader/ACF.dat')
        self.bader_atoms = deepcopy(self.atoms)

        if acf.is_file():
            self.bader = True
            attach_charges(self.bader_atoms, self.homedir + '/ACF.dat')
        elif acf_folder.is_file():
            self.bader = True
            attach_charges(self.bader_atoms, self.homedir + '/bader/ACF.dat')



