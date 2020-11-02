
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
    specifically for adsorbates 

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

    def __init__(self, findiffdb, fields_to_choose,\
            direction, displacement, atomsIS, atomsFS, db_params={}):

        self.findiffdb = findiffdb # database of results with finite diff calcualtions
        self.fields_to_choose = fields_to_choose
        self.db_params = db_params

        self.vibresults = {} # Results from the vibrational calculation 
        self.dmudR = {} # change of dipole moment accompanying disp of atoms 
        self.dFdG = {} # Change in the force with the electric field applied
        self.direction = direction
        self.displacement = displacement
        self.atomsIS = atomsIS
        self.atomsFS = atomsFS
        self.indices = []
        self.modes = []
        self.q = {}

        # Default get the information from the forces and the electic field 
        self._data_for_force()
        # self._get_reaction_path()



        
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
            indice = int(indice_str.split(direction)[0])
            try:
                disp = float(row.displacement.split('_')[-2])
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
            self.vibresults.setdefault(state,{})\
                .setdefault(indice,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['energy'] = row.energy
            self.vibresults.setdefault(state,{})\
                .setdefault(indice,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['forces'] = forces
            self.vibresults.setdefault(state,{})\
                .setdefault(indice,{}).setdefault(disp,{}).setdefault(direction,{}).setdefault(field,{})['dipole'] = dipole
    
    
    def get_dmudR(self):
        for state in self.vibresults:
            for indice in self.vibresults[state]:

                try:
                    mu_p = self.vibresults[state][indice][self.displacement]['p'][0.0]['dipole']
                    mu_m = self.vibresults[state][indice][self.displacement]['m'][0.0]['dipole']
                except KeyError:
                    continue

                delta_mu = mu_p - mu_m
                dmudR = delta_mu / 2 / self.displacement
                self.dmudR.setdefault(state,{})[indice] = dmudR


    def get_dFdG(self):

        for state in self.vibresults:
            for indice in self.vibresults[state]:

                self.indices.append(indice)

                fields_forces = {}
                fields = []

                try:
                    for field in self.vibresults[state][indice][self.displacement][self.direction]:
                        fields_forces[field] = self.vibresults[state][indice][self.displacement][self.direction][field]['forces']
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


    def get_reaction_path(self):
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