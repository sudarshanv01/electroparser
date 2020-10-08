
import numpy as np
from ase import units
import ase
from pprint import pprint
from matplotlib import cm
from useful_functions import get_fit_from_points
import matplotlib.pyplot as plt
from ase.utils import pickleload
import os

class ForceExtrapolationData:

    """
    Extrapolation based on the forces at a finite difference electric fields 

    ...

    Attributes
    ----------

    findiffdb : ase db 
        Finite difference database 
    fields_to_choose : list 
        Fields over which the finite differencing is done
    db_params : dict 
        Parameters which will be fed into the ase db select option

    vibresults : dict 
        Results stored from finite difference calculations 
    dmudR : dict 
        Change in dipole moment with finite diff dispacement 
    dFdG : dict 
        Change in the Forces along with application of Electric field
    
    Methods
    -------

    get_dmudR
        Gets the change in the dipole moment along z directions ( periodic calcs )
    get_dFdG 
        Gets the change in the forces along application of the electric field
    

    """

    def __init__(self, findiffdb, fields_to_choose, db_params={}):

        self.findiffdb = findiffdb # database of results with finite diff calcualtions
        self.fields_to_choose = fields_to_choose
        self.db_params = db_params

        self.vibresults = {} # Results from the vibrational calculation 
        self.dmudR = {} # change of dipole moment accompanying disp of atoms 
        self.dFdG = {} # Change in the force with the electric field applied


        # Default get the information from the forces and the electic field 
        self._data_for_force()

        
    def _data_for_force(self):
        # get the data for extrapolating through dipoles and or forces 
        for row in self.findiffdb.select(**self.db_params):
            try:
                conc = row.proton_conc
            except AttributeError:
                conc = 'conc'
            try:
                cell = row.cell_size
            except AttributeError:
                cell = 'cell'
            try:
                state = row.states.replace('_SP', '')
            except AttributeError:
                state = 'state'
            indice_str = row.findiff # The indice that was moved and in which direction 
            direction = 'p' if 'p' in indice_str else 'm'
            field = row.field
            indice = int(indice_str.split(direction)[0])
            try:
                disp = float(row.sampling.split('_')[-2])
            except ValueError:
                continue
            atoms = row.toatoms()
            try:
                # this is where SJM stores the dipole moment
                dipole = atoms.get_dipole_moment()
            except ase.calculators.calculator.PropertyNotImplementedError:
                # Not an sjm calculation
                dipole = row.dipole_field
            forces = atoms.get_forces()[indice]
            self.vibresults.setdefault(conc,{}).setdefault(cell,{}).setdefault(state,{})\
                .setdefault(indice,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['energy'] = row.energy
            self.vibresults.setdefault(conc,{}).setdefault(cell,{}).setdefault(state,{})\
                .setdefault(indice,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['forces'] = forces
            self.vibresults.setdefault(conc,{}).setdefault(cell,{}).setdefault(state,{})\
                .setdefault(indice,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['dipole'] = dipole
    
    
    def get_dmudR(self):
        for conc in self.vibresults:
            for cell in self.vibresults[conc]:
                for state in self.vibresults[conc][cell]:
                    for indice in self.vibresults[conc][cell][state]:
                        
                        for disp in  self.vibresults[conc][cell][state][indice]:

                            # Take the dipole moment difference 
                            # For the case of zero field 

                            try:
                                mu_p = self.vibresults[conc][cell][state][indice][disp]['p'][0.0]['dipole']
                                mu_m = self.vibresults[conc][cell][state][indice][disp]['m'][0.0]['dipole']
                            except KeyError:
                                continue

                            delta_mu = mu_p - mu_m
                            dmudR = delta_mu / 2 / disp
                            self.dmudR.setdefault(conc,{}).setdefault(cell,{}).setdefault(state,{})\
                                .setdefault(disp,{})[indice] = dmudR


        
    
    def get_dFdG(self):

        for conc in self.vibresults:
            for cell in self.vibresults[conc]:
                for state in self.vibresults[conc][cell]:
                    for indice in self.vibresults[conc][cell][state]:
                        
                        for disp in  self.vibresults[conc][cell][state][indice]:
                            for direc in self.vibresults[conc][cell][state][indice][disp]:

                                fields_forces = {}
                                fields = []

                                # pprint(self.vibresults[conc])

                                for field in  self.vibresults[conc][cell][state][indice][disp][direc]:
                                    fields_forces[field] = self.vibresults[conc][cell][state][indice][disp][direc][field]['forces']
                                    fields.append(field)
                                if not self.fields_to_choose:
                                    f_pos = sorted([ff for ff in fields if ff > 0.0]) 
                                else:
                                    # fields are provided
                                    f_pos = self.fields_to_choose
                                fmax = np.max(f_pos) ; fmin = np.min(f_pos)
                                deltaf = fmax - fmin
                                dFdG = ( -1 * fields_forces[fmax] + 8 * fields_forces[fmin] \
                                            - 8 * fields_forces[-1*fmin] + fields_forces[-1*fmax] ) / 6.0 / 2 / (deltaf)

                                try:
                                    self.dFdG.setdefault(conc,{}).setdefault(cell,{})\
                                        .setdefault(state,{}).setdefault(direc,{})[disp].append(dFdG)
                                except (AttributeError, KeyError):
                                    self.dFdG.setdefault(conc,{}).setdefault(cell,{})\
                                        .setdefault(state,{}).setdefault(direc,{})[disp] = []
                                    self.dFdG.setdefault(conc,{}).setdefault(cell,{})\
                                        .setdefault(state,{}).setdefault(direc,{})[disp].append(dFdG)

                                    
class EigenModesHessian:
    """
    
    Collects the Hessian and Eignemodes for a given 
    transition state calculation

    :param directory: Directory where the pickle files are stored
    :type directory: str
    :atoms: Atoms object which corresponds to the transition state 
    :type atoms: atoms object

    """

    def __init__(self, atoms, directory, indices, dx, prefix):
        self.atoms = atoms # atoms object for TS
        self.directory = directory # directory where the calculation was done 
        self.indices = indices # indices making up the vibrations modes
        self.dx = dx # displacement that the atoms had
        self.prefix = prefix # prefix assigned to vibration file

        self.H = np.empty((3*len(self.indices),3*len(self.indices))) # Dynamical matrix

        self.modes = [] # eigenmodes 
        self.frequencies = [] # eigenfrequencies
        self.frequencies_cm = [] # Real frequencies in cm-1
        self.unit_vectors = [] # Unit vectors to perturb atoms along normal mode

        # get the Hessian
        self.Hessian()
        self.eigenmodes()

    def Hessian(self):
        # Get the equilibrium forces 
        forces0 = pickleload(open(os.path.join(self.directory,self.prefix+'.eq.pckl'), 'rb'))
        # Gets the Hessian from the pickle file
        r = 0
        for i, index in enumerate(self.indices):
            for axes in 'xyz': # Moved along all direction 
                pickle_filename = f'{self.prefix}.{index}{axes}'
                forces_p = pickleload(open(os.path.join(self.directory,pickle_filename+'+.pckl'), 'rb'))
                forces_n = pickleload(open(os.path.join(self.directory,pickle_filename+'-.pckl'), 'rb'))
                self.H[r] = (forces_n - forces_p)[self.indices].ravel() / 4. / self.dx
                r += 1
    

    def eigenmodes(self):
        self.H += self.H.copy().T
        m = self.atoms.get_masses()[self.indices]
        self.im = np.repeat(m**-0.5, 3)
        omega2, modes = np.linalg.eigh(self.im[:,None] *  self.H * self.im )
        self.modes = modes
        s = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
        self.frequencies = s * omega2.astype(complex)**0.5
        self.frequencies_cm = np.real(0.01 * units._e / units._c / units._hplanck * self.frequencies)
