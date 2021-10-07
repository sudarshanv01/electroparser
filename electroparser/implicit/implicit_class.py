#!/usr/bin/python


""" Class to parse and plot implicit solvent calculations """

from useful_functions import AutoVivification, get_vasp_nelect0
import os
from ase.data import atomic_numbers
import numpy as np


class ImplicitParse:

    """ 
    Reads in the ASE database and gives all the information 
    needed to do electrochemistry 
    if an implicit solvent is used 
    """

    def __init__(self, electrodb, config):
        self.electrodb = electrodb # Database to parse
        self.config = config # Config to query database

        self.structures = AutoVivification() # All ase structures
        self.dft_workfunctions = AutoVivification() # workfunction for all calcs
        self.other_info = AutoVivification()

        self.selected_rows = []

        self.get_information()


    # Standard querying of database
    def _query_electrodb(self, kwargs):
         # Function to query an ase database with the needed kwargs
         for row in self.electrodb.select(**kwargs):
             self.selected_rows.append(row)

    def _get_surface_area(self, ase_atoms):
        height = ase_atoms.get_cell()[-1,-1]
        volume = ase_atoms.get_volume()
        area = volume / height
        return area


    def get_information(self):
         config = self.config
         self._query_electrodb(config)

         for row in self.selected_rows:

             ase_atoms = row.toatoms()
             if row.paw == 'vasp':
                 nelect0 = get_vasp_nelect0(ase_atoms)
                 nelect = row.tot_charge
                 charge = round(nelect - nelect0, 2)
             else:
                 # This is a QE run
                 charge = round(row.charge, 2)

             # Must have state that is non negotiable
             state = row.states.replace('IS_SP', 'IS').replace('FS_SP', 'FS')
             try:
                 # If there is a proton
                 chosen_proton = row.proton_conc
             except AttributeError:
                 chosen_proton = 'nan'
             try:
                 cell = row.cell_size
             except:
                 cell = 'nan'

             surface_area = row.volume / row.cell[-1,-1]


             self.dft_workfunctions[chosen_proton][cell][state][charge] = row.wf
             self.other_info[chosen_proton][cell][state][charge]['surface_area'] = self._get_surface_area(ase_atoms)

             try:
                 # Get the energy from the row
                 self.other_info[chosen_proton][cell][state][charge]['energy'] = row.energy
             except AttributeError:
                 # Get the bader charge from the file
                 atoms = row.toatoms()
                 # print([atom.symbol for atom in atoms])
                 #print([atomic_numbers[atom.symbol] - atom.charge for atom in atoms])
                 count_index = []
                 if row.states == 'FS':
                     # exlude the adsobed H
                     Pt_indices = [atom.index for atom in atoms if atom.symbol == 'Pt']
                     OH_indices = [atom.index for atom in atoms if atom.symbol in ['O', 'H']]
                     adsorbed_H = np.argpartition([ np.min(atoms.get_distances(a, Pt_indices)) for a in OH_indices ],2)[0:2]
                     OH_indices = np.delete(OH_indices, adsorbed_H)
                 else:
                    OH_indices = [atom.index for atom in atoms if atom.symbol in ['O', 'H']]
                 # Get the total charge
                 charge_paw = {'Pt':10, 'H':1, 'O':6}
                 total_charge = np.sum([ [charge_paw[atom.symbol] - atomic_numbers[atom.symbol] - atom.charge for atom in atoms[OH_indices]] ])
                 self.other_info[chosen_proton][cell][state][charge]['bader'] = total_charge / 2

class FracChargeParse:

    def __init__(self, electrodb, config):
        self.electrodb = electrodb # Database to parse
        self.config = config # Config to query database

        self.structures = AutoVivification() # All ase structures
        self.dft_workfunctions = AutoVivification() # workfunction for all calcs
        self.other_info = AutoVivification()

        self.selected_rows = []

        self.get_information()


    # Standard querying of database
    def _query_electrodb(self, kwargs):
         # Function to query an ase database with the needed kwargs
         for row in self.electrodb.select(**kwargs):
             self.selected_rows.append(row)

    def _write_out_atoms(self, atoms, name, folder):
        os.system('mkdir -p ' + folder)
        atoms.write(folder + '/' + name + '.traj')

    def _get_surface_area(self, ase_atoms):
        height = ase_atoms.get_cell()[-1,-1]
        volume = ase_atoms.get_volume()
        area = volume / height
        return area


    def get_information(self):
         config = self.config
         self._query_electrodb(config)
         for row in self.selected_rows:

             ase_atoms = row.toatoms()
             nelect0 = get_vasp_nelect0(ase_atoms)
             nelect = row.tot_charge
             if row.paw == 'vasp':
                 try:
                     charge = float(row.structure.split('_')[-1].replace('p', '.'))
                 except AttributeError:
                     continue
             else:
                 # This is a QE run
                 charge = round(nelect, 2)

             # Must have state that is non negotiable
             state = row.states.replace('IS_SP', 'IS').replace('FS_SP', 'FS')
             try:
                 # If there is a proton
                 chosen_proton = row.proton_conc
             except AttributeError:
                 chosen_proton = 'nan'
             try:
                 cell = row.cell_size
             except:
                 cell = 'nan'

             surface_area = row.volume / row.cell[-1,-1]

             # Write out the atoms
             self._write_out_atoms(ase_atoms, str(cell) + '_' + state + '_charge_' + str(round(charge,2)) + '_conc_' + chosen_proton, 'output_structures/implicit')

             self.dft_workfunctions[chosen_proton][cell][state][charge] = row.wf
             self.other_info[chosen_proton][cell][state][charge]['surface_area'] = self._get_surface_area(ase_atoms)
             self.other_info[chosen_proton][cell][state][charge]['energy'] = row.energy
