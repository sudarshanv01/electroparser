#!/usr/bin/python


""" Helper functions to parse """

import os
import ase
from copy import deepcopy
import numpy as np
from pprint import pprint

# Gets the needed directories directly under the homedir given
def _get_dirs(homedir):
    # Given a path, returns a list of directories
    # Also returns if state information is provided

    directories = []
    dont_look_here = ['noparse', 'cineb', 'neb', 'archive', 'initial_structures', 'template_scripts', 'output', 
            ]

    for direct in os.listdir(homedir):
        if os.path.isdir(homedir + '/' + direct) and direct not in dont_look_here:
            if 'cineb' not in direct and 'neb' not in direct:
                directories.append(homedir + '/' +  str(direct) )
    return directories

# Get the lowest directory
# TODO: This should be one line
def _get_homedir(levels):
    # INPUTS
    # levels to parse from

    paths = {}

    for i in range(levels):
        level = i
        if level == 0:
            current_list = _get_dirs('.')
            paths[level] = current_list

        elif level > 0:
            accumulated_list = []
            for j in range(len(paths[level-1])):
                current_list = _get_dirs(paths[level-1][j])
                accumulated_list.append(current_list)
            saved_list = [item for sublist in accumulated_list for item in sublist]
            paths[level] = saved_list

    print("Paths explored:")
    pprint(paths)
    homedirs = paths[levels-1]

    return homedirs

# Hard coded set of instructions based on the my directory structure
def _initialize_specifics(level):
    # Get the homedir
    homedirs = _get_homedir(level)

    # Get states
    states = []
    facet = []

    # possibities to check when parsing a directory
    states_possible = ['state_', 'states_', 'slab','CO_site', 'CH', 'Li', 'Pb', 'IS', 'FS', \
            'wf', 'neb', 'cineb', 'TS', 'OH', 'CO2']
    facets_possible = ['DV', 'SW', 'SV', 'DV4N', 'facet']
    cell_size_possible = ['x']
    structure_possible = ['structure', 'image_', 'theta_']
    sampling_possible = ['sampling', 'ps', 'disp_']
    termination_possible = ['termination']
    proton_conc_possible = ['_conc']
    dopant_p_block_possible = ['_doped']
    vacancy_number_possible = ['vacancy_']
    dopant_number_possible = ['dopant_']
    findiff_possible = ['pz', 'mz', 'py', 'my', 'px', 'mx']
    metal_dopant_possible = ['metal_dopant']
    possibilities = {'states':states_possible,
                     'facets':facets_possible,
                     'cell_size':cell_size_possible,
                     'sampling':sampling_possible,
                     'termination':termination_possible,
                     'proton_conc':proton_conc_possible,
                     'structure':structure_possible,
                     'termination':termination_possible,
                     'dopant_p_block':dopant_p_block_possible,
                     'vacancy_number':vacancy_number_possible,
                     'dopant_number':dopant_number_possible,
                     'findiff':findiff_possible,
                     'metal_dopant':metal_dopant_possible,
                    }
    states = np.nan
    facet = np.nan
    cell_size_input = np.nan
    structure_input = np.nan
    sampling_type = np.nan
    termination_type = np.nan
    species_type = np.nan
    proton_conc_type = np.nan
    sampling_type = np.nan

    variable_poss = {'states':states,
                     'facets':facet,
                     'cell_size':cell_size_input,
                     'structure':structure_input,
                     'termination':termination_type,
                     'proton_conc':proton_conc_type,
                     'sampling':sampling_type,
                    }
    return {'possibilities':possibilities, 'variable_poss':variable_poss,
                'homedirs':homedirs}

def check_if_complete():
    # reads in .complete_data
    completed_files = []
    try:
        with open('.complete_data', 'r') as f:
            lines = f.readlines()
            for line in lines:
                completed_files.append(line.split('/\n')[0])
    except FileNotFoundError:
        completed_files = []

    return completed_files


def data_to_store(level, db):
    from parser_class import Parser
    from neb_parser_class import NebParser
    initialize = _initialize_specifics(level)
    possibilities = initialize['possibilities']
    variable_poss = initialize['variable_poss']
    homedirs_all = initialize['homedirs']
    # Check if a previous parsing has occured
    already_complete = check_if_complete()
    homedirs = [ fil for fil in homedirs_all if fil not in already_complete ]


    for i in range(len(homedirs)):
        data = {}

        homedir = homedirs[i] + '/'
        # Get information directly from the directories
        split_list = homedir.split('/')
        for possibility in possibilities:
            for check in possibilities[possibility]:
                for lst in split_list:
                    if check in lst :
                        data[possibility] = lst
        pprint(data)
        # Check if parsing a neb or a static calculation
        try:
            if 'neb' in  data['states']:
                print('Parsing neb from folder:' + homedir)
                neb_run = True
                base = data['states']
            else:
                print('Parsing folder:'  + homedir)
                neb_run = False
        except KeyError:
            neb_run = False

        # Invoke and instance of the Parser class
        if neb_run:
            parse = NebParser(homedir)
            metadata = {}

            atoms = parse.atoms
            if not bool(atoms):
                # No atoms object in this directory
                print('Skipping')
            else:
                for index, atom in enumerate(atoms):
                    data['atoms'] = atom
                    data['wf'] = parse.wf[index]
                    data['implicit'] = parse.implicit
                    data['field'] = parse.field
                    data['dipole_field'] = parse.dipole[index]
                    data['tot_charge'] = parse.charge
                    metadata['input_file'] = parse.input_file
                    data['paw'] = parse.paw
                    data['functional'] = parse.functional
                    metadata['ldau'] = parse.ldau_data
                    data['states'] = base + str(index)
                    #data['ir_intensity'] = parse.ir_intensity
                    #data['sigma_raman'] = parse.raman_intensity

                    db.write(**data)


        else:
            try:
                parse = Parser(homedir)
            except IndexError:
                print('Skipping '  + homedir)
                continue
            except ValueError:
                print('Skipping {s}, something wrong with the WF parsing'.format(s=homedir))
                continue
            except AssertionError:
                print('Skipping {s}, something wrong with the QE output'.format(s=homedir))
                continue
            except ase.io.formats.UnknownFileTypeError:
                print('Error reading ASE atoms file')
                continue
            metadata = {}

            atoms = parse.atoms
            if not bool(atoms):
                # No atoms object in this directory
                print('Skipping')
            else:
                data['atoms'] = atoms
                data['wf'] = parse.wf
                data['implicit'] = parse.implicit
                data['field'] = parse.field
                data['pw'] = parse.pw
                data['dipole_field'] = parse.dipole
                data['tot_charge'] = parse.charge
                metadata['input_file'] = parse.input_file
                data['paw'] = parse.paw
                data['functional'] = parse.functional
                data['ldau'] = parse.ldau
                metadata['vibrations'] = parse.vibrations
                if data['ldau']:
                    # This is a LDA+U calculation 
                    metadata['ldau'] = parse.ldau_data
                # Only save d-band data for slab calculations
                if isinstance(parse.spin_pol, bool) and 'slab' in data['states']:
                    if parse.spin_pol:
                        data['d_centre_up'] = parse.d_centre_up
                        data['d_centre_down'] = parse.d_centre_down
                        data['max_hilbert_up'] = parse.max_hilbert_up
                        data['max_hilbert_down'] = parse.max_hilbert_down
                        data['occupancy_up'] = parse.occupancy_up
                        data['occupancy_down'] = parse.occupancy_down
                    else:
                        data['d_centre'] = parse.d_centre
                        data['max_hilbert'] = parse.max_hilbert
                        data['occupancy'] = parse.occupancy

                # Update energy if only POSCARS stored
                if bool(parse.energy):
                    data['energy_parsed'] = parse.energy
                data['data'] = metadata

                db.write(**data)

                if parse.bader:
                    # There is a bader atoms object 
                    data_bader = deepcopy(data)
                    data_bader['atoms'] = parse.bader_atoms
                    db.write(**data_bader)

        # On correct parsing write out to .complete_data
        with open('.complete_data', 'a') as f:
            f.write(homedir)
            f.write('\n')
