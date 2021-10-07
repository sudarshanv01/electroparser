
import numpy as np
from ase import units
import ase
from pprint import pprint
from matplotlib import cm
from useful_functions import get_fit_from_points
import matplotlib.pyplot as plt
from ase.utils import pickleload
import os
from ase.geometry import geometry
from pprint import pprint
class ForceExtrapolation:

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

    def __init__(self, findiffdb, fields_to_choose, indices, atomsIS, atomsFS,\
            direction, displacement, conc, cell_size, db_params={}):

        self.findiffdb = findiffdb # database of results with finite diff calcualtions
        self.fields_to_choose = fields_to_choose
        self.db_params = db_params
        self.atomsIS = atomsIS
        self.atomsFS = atomsFS
        self.indices = indices

        self.vibresults = {} # Results from the vibrational calculation 
        self.dmudR = {} # change of dipole moment accompanying disp of atoms 
        self.dFdG = {} # Change in the force with the electric field applied
        self.q = {} # Get the charge based on dFdG
        self.modes = [] # Modes for doing the dot product
        self.direction = direction
        self.displacement = displacement
        self.conc = conc
        self.cell_size = cell_size

        if bool(self.conc):
            self.db_params['proton_conc'] = self.conc
        self.db_params['cell_size'] =  self.cell_size

        # Default get the information from the forces and the electic field 
        self._data_for_force()
        self._get_reaction_path()



        
    def _data_for_force(self):
        # get the data for extrapolating through dipoles and or forces 
        for row in self.findiffdb.select(**self.db_params):
            # Populate the results dictionary based on some presets
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
            indice = [int(a) for a in indice_str[:-2].split('_')] #int(indice_str.split(direction)[0])
            indice_s = indice_str[:-2]

            try:
                disp = float(row.sampling.split('_')[-2])
            except ValueError:
                continue
            except AttributeError:
                disp = float(row.displacement.split('_')[-2])

            atoms = row.toatoms()

            try:
                # this is where SJM stores the dipole moment
                dipole = atoms.get_dipole_moment()
            except ase.calculators.calculator.PropertyNotImplementedError:
                # Not an sjm calculation
                dipole = row.dipole_field

            forces = atoms.get_forces()[indice]
            self.vibresults.setdefault(state,{})\
                .setdefault(indice_s,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['energy'] = row.energy
            self.vibresults.setdefault(state,{})\
                .setdefault(indice_s,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['forces'] = forces
            self.vibresults.setdefault(state,{})\
                .setdefault(indice_s,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['dipole'] = dipole
    
    
    def get_dmudR(self):
        for indice in self.indices:
            for state in self.vibresults:

                # try:
                mu_p = self.vibresults[state][str(indice)][self.displacement]['p'][0.0]['dipole']
                mu_m = self.vibresults[state][str(indice)][self.displacement]['m'][0.0]['dipole']
                # except KeyError:
                #     continue

                delta_mu = mu_p - mu_m
                dmudR = delta_mu / 2 / self.displacement
                self.dmudR.setdefault(state,{})[indice] = dmudR


    def get_dFdG(self):

        for indice in self.indices:
            for state in self.vibresults:

                fields_forces = {}
                fields = []

                try:
                    for field in self.vibresults[state][str(indice)][self.displacement][self.direction]:
                        fields_forces[field] = self.vibresults[state][str(indice)][self.displacement][self.direction][field]['forces']
                        fields.append(field)
                except KeyError:
                    continue

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
                    self.dFdG[state].append(dFdG)
                except (AttributeError, KeyError):
                    self.dFdG[state] = []
                    self.dFdG[state].append(dFdG)

    def get_q(self):
        """ Get the q by taking the dot product between the dFdG and the 
        reaction mode 
        """
        for state in self.vibresults:
            dFdG = []
            j = 0
            for i in range(3*len(self.indices)):
                if (i+1)%3 == 0:
                    # a z-component
                    try:
                        differential = self.dFdG[state][j]
                    except IndexError:
                        print('Missing data!')
                        continue
                    dFdG.append([0,0,differential[-1]])
                    j += 1
                else:
                    dFdG.append([0, 0, 0])
            dFdG = np.array(dFdG)
            mu_axes = dFdG.T[-1]
            # now dot product with the different modes available
            for index, mode in enumerate(self.modes):
                try:
                    q = np.dot(mu_axes, mode)
                except ValueError:
                    continue
                self.q.setdefault(state,{})[index] = q




    def _get_reaction_path(self):
        """Gets the reaction path from the atoms object
        """
        ## check if the atoms are on the same side of the unit cell
        cell = self.atomsIS.get_cell() # same cell used in IS and FS hopefully
        # get the vector respresenting the difference of the two 
        vector_all = self.atomsIS.get_positions() - self.atomsFS.get_positions()
        vectors = vector_all[self.indices]
        min_vec = []
        for v in vectors:
            vmin, vlen = geometry.find_mic(v, cell, pbc=True)
            min_vec.append(vmin)
        ravel_vec = np.ravel(min_vec)
        self.modes.append( ravel_vec / np.linalg.norm(ravel_vec) )
    


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
