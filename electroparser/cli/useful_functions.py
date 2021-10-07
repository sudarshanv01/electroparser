#!/usr/bin/python

from ase.io import read, write

def get_reference_energies():
    # results a dictionary with all the gas phase energies 
    # static link to database
    results = {}
    from ase.db import connect 
    database_file = '/home/cat/vijays/project/gas_phase_energies/gas_phase.db'
    database = connect(database_file)
    for row in database.select():
        state = row.states
        pw = row.pw 
        functional = row.functional
        results.setdefault(state,{}).setdefault(functional,{}).setdefault(pw,{})['energy'] = row.energy
        results[state][functional][pw]['vibrations'] = row.data.vibrations
        results[state][functional][pw]['atoms'] = row.toatoms()
    
    return results

        
def get_plot_params():
    """Create the plot parameters used in the plotting
    all the figures in the paper
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.use('Agg')

    mpl.rcParams['font.family'] = 'DejaVu Serif'
    # mpl.rcParams['font.serif'] = 'Times New Roman'
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams.update({
    # "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"]})
    })

    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 2
    plt.rcParams['xtick.minor.size'] = 7
    plt.rcParams['xtick.minor.width'] = 2
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 2
    plt.rcParams['ytick.minor.size'] = 7
    plt.rcParams['ytick.minor.width'] = 2

    plt.rcParams['axes.labelsize'] = 28

def get_vasp_nelect0(atoms):
    import pickle
    import numpy as np
    import os
    # Given an atoms object it will return the nelect at zero charge
    # Based on default vasp potentials
    with open('/p/home/jusers/vijay1/juwels/bin/scripts/nelect0.pickle','rb') as handle:
        nelect0 = pickle.load(handle)
    default_nelect = []
    for i in range(len(atoms)):
        default_nelect.append(nelect0[atoms[i].symbol])
    default_nelect = np.array(default_nelect)
    return default_nelect.sum()


def get_fit_from_points(x, y, order):
    import numpy as np
    fit = np.polyfit(x, y, order)
    p = np.poly1d(fit)
    return {'fit':fit, 'p':p}


def read_dos(homedir):
    import pickle
    f = open(homedir + '/dos.pickle', 'rb')
    energies, dos, pdos = pickle.load(f, encoding='latin1')
    f.close()
    return energies, dos, pdos

# Correct for energies
# Espresso takes in only the Fermi energy from the old log file
# Need to change this in the future
def correction(homedir):
    import os
    old_log = os.popen('grep -a Fermi '+ homedir + '/log |tail -2', 'r')
    efermi_old = float(old_log.readline().split()[-2])
    new_log = os.popen('grep -a Fermi ' + homedir + '/pwnscf.log |tail -1', 'r')
    efermi_new = float(new_log.readline().split()[-2])

    return efermi_new - efermi_old

def mean_d_band(energies, pdos):
    from scipy.integrate import quad, simps
    def integrand_numerator(e, rho):
        return e * rho
    def integrand_denominator(e,rho):
        return rho
    epsilon_d_num = simps(integrand_numerator(energies,pdos), energies)
    epsilon_d_den = simps(integrand_denominator(energies,pdos), energies)

    return epsilon_d_num / epsilon_d_den

def new_fermi_level(homedir):
    import os
    new_log = os.popen('grep -a Fermi ' + homedir + '/pwnscf.log |tail -1', 'r')
    efermi_new = float(new_log.readline().split()[-2])
    return efermi_new


def surface_area(homedir):
    from ase.io import read
    atoms = read(homedir + '/spe.traj')
    height = atoms.get_cell()[1, 1]
    volume = atoms.get_volume()
    sa = volume / height
    return sa


def quantum_capacitance_FBA(sa, energies, fermi, dos):
    import numpy as np
    # Takes list of energies, fermi and dos and gives out
    # the Quantum capacitance based on Wood (2015)
    from scipy.integrate import simps
    from ase.units import Rydberg
    e = 1.6e-19 # C
    micro = 1e6
    A2tocm2 = 1e-16
    # Quantum capacitance is e^2 n(mu,N)
    closest_to_zero = min(abs(energies))
    position_fermi = np.where(energies == closest_to_zero)[0]
    if position_fermi.size == 0:
        position_fermi = np.where(energies == -1 * closest_to_zero)[0]

    Cq = ( dos[0][position_fermi] + dos[1][position_fermi] ) * ( 1 / sa ) * ( e * micro / A2tocm2 )
    return Cq

def get_bands_DOS(pdos):
    if len(pdos) == 10: # spin unpolarised
        s, p_y, p_z, p_x, d_xy, d_yz, d_z2, d_xz, d_x2_y2 =  pdos[1:,1:]
        data = {'s':s, 'p_y':p_y, 'p_z':p_z, 'p_x':p_x, 'd_xy':d_xy, 'd_yz':d_yz,\
                'd_z2':d_z2, 'd_xz':d_xz, 'd_x2_y2':d_x2_y2}
        return data

    elif len(pdos) == 19: # spin polarized
        s_up, s_down, p_y_up, p_y_down, p_z_up, p_z_down, p_x_up, p_x_down, \
                d_xy_up, d_xy_down, d_yz_up, d_yz_down, d_z2_up, d_z2_down, \
                d_xz_up, d_xz_down, d_x2_y2_up, d_x2_y2_down =  pdos[1:,1:]
        data = {'s_up':s_up, 's_down':s_down, 'p_y_up':p_y_up, 'p_y_down':p_y_down,\
                'p_z_up':p_z_up, 'p_z_down':p_z_down, 'p_x_up':p_x_up, 'p_x_down':p_x_down,\
                'd_xy_up':d_xy_up, 'd_xy_down':d_xy_down, 'd_yz_up':d_yz_up, 'd_yz_down':d_yz_down, \
                'd_z2_up':d_z2_up, 'd_z2_down':d_z2_down, 'd_x2_y2_up':d_x2_y2_up, 'd_x2_y2_down':d_x2_y2_down,\
                'd_xz_up':d_xz_up, 'd_xz_down':d_xz_down}
        return data
    print("Missing data from DOSCAR - check calculation")


def d_band_info(energies, pdos):
    from scipy.integrate import quad, simps
    def integrand_numerator(e, rho):
        return e * rho
    def integrand_denominator(e,rho):
        return rho
    # Getting mean d-band
    epsilon_d_num = simps( integrand_numerator(energies,pdos), energies)
    epsilon_d_den = simps(integrand_denominator(energies,pdos), energies)
    epsilon_d = epsilon_d_num / epsilon_d_den
    # First moment of d-band
    def integrand_moment_numerator(e, epsilon_d, moment, rho):
        return ( (e - epsilon_d ) ** moment ) * rho
    w_d_numerator = simps( integrand_moment_numerator(energies, epsilon_d, 2, \
            pdos), energies)
    w_d_denominator = epsilon_d_num
    w_d = w_d_numerator / w_d_denominator

    return [epsilon_d, w_d]


def get_filled_states(energies, pdos):
    # does the integral of -inf to Ef for the pdos
    from scipy.integrate import quad, simps
    energy_fill = []
    pdos_occ = []
    for index, energy in enumerate(energies):
        if energy < 0:
            energy_fill.append(energy)
            pdos_occ.append(pdos[index])
    def integrand(e, rho):
        return rho

    f = simps(integrand(energy_fill, pdos_occ), energy_fill)

    return f





def adsorbates(atoms, ads_name, position, h):
    from ase.build import molecule, add_adsorbate
    new_atoms = atoms.copy()
    if ads_name == 'COa':
        co = molecule('CO')
        add_adsorbate(new_atoms, co, h, position=position[0:2], mol_index=1)
    elif ads_name == 'COg':
        co = molecule('CO')
        add_adsorbate(new_atoms, co, h, position=position[0:2], mol_index=1)
    elif ads_name == 'COOHa':
        hcooh = molecule('HCOOH')
        cooh = hcooh.copy()
        del cooh[-1]
        cooh.rotate(90,'x')
        cooh.rotate(180,'y')
        cooh.rotate(180, 'z')
        add_adsorbate(new_atoms, cooh, h, position=position[0:2], mol_index=1)
    elif ads_name == 'CO2a':
        co2 = read('/home/cat/vijays/tools/adsorbates/co2.xyz')
        add_adsorbate(new_atoms, co2, h, position=position[0:2], mol_index=0)
    elif ads_name == 'CO2g':
        co2 = molecule('CO2')
        co2.rotate(90, 'x')
        add_adsorbate(new_atoms, co2, h, position=position[0:2], mol_index=0)
    elif ads_name == 'H':
        H = Atoms('H')
        add_adsorbate(new_atoms, H, h, position=position[0:2])
    return new_atoms


def read_catmap_table(txt):
    f = open(txt)
    lines = f.readlines()
    values = []
    header = lines[0].split()
    print(header)
    del lines[0]
    for line in lines:
        values.append(line.split())
    values_dict = {}
    for i in range(len(header)):
        values_dict[header[i]] = [ float(value[i]) for value in values ]

    return values_dict


def add_adsorbates(atoms, ads_name, position, h):
    from ase.build import molecule, add_adsorbate
    from ase import Atoms
    new_atoms = atoms.copy()
    if ads_name == 'CO':
        co = molecule('CO')
        add_adsorbate(new_atoms, co, h, position=position[0:2], mol_index=1)
    elif ads_name == 'COg':
        co = molecule('CO')
        add_adsorbate(new_atoms, co, h, position=position[0:2], mol_index=1)
    elif ads_name == 'COOH':
        hcooh = molecule('HCOOH')
        cooh = hcooh.copy()
        del cooh[-1]
        cooh.rotate(90,'x')
        cooh.rotate(180,'y')
        cooh.rotate(180, 'z')
        add_adsorbate(new_atoms, cooh, h, position=position[0:2], mol_index=1)
    elif ads_name == 'CO2':
        co2 = read('/home/cat/vijays/tools/adsorbates/co2.xyz')
        add_adsorbate(new_atoms, co2, h, position=position[0:2], mol_index=0)
    elif ads_name == 'CHO':
        cho = read('/home/cat/vijays/tools/adsorbates/CHO.traj')
        cind = [ atom.index for atom in cho if atom.symbol == 'C']
        add_adsorbate(new_atoms, cho, h, position=position[0:2], mol_index=cind[0])
    elif ads_name == 'CH2O':
        ch2o = read('/home/cat/vijays/tools/adsorbates/CH2O.traj')
        cind = [ atom.index for atom in ch2o if atom.symbol == 'C']
        add_adsorbate(new_atoms, ch2o, h, position=position[0:2], mol_index=cind[0])
    elif ads_name == 'CH2OH':
        ch2oh = read('/home/cat/vijays/tools/adsorbates/CH2OH.traj')
        cind = [ atom.index for atom in ch2oh if atom.symbol == 'C']
        add_adsorbate(new_atoms, ch2oh, h, position=position[0:2], mol_index=cind[0])
    elif ads_name == 'CH2':
        ch2 = read('/home/cat/vijays/tools/adsorbates/CH2.traj')
        cind = [ atom.index for atom in ch2 if atom.symbol == 'C']
        add_adsorbate(new_atoms, ch2, h, position=position[0:2], mol_index=cind[0])
    elif ads_name == 'CH3':
        ch3 = read('/home/cat/vijays/tools/adsorbates/CH3.traj')
        cind = [ atom.index for atom in ch3 if atom.symbol == 'C']
        add_adsorbate(new_atoms, ch3, h, position=position[0:2], mol_index=cind[0])
    elif ads_name == 'CO2g':
        co2 = molecule('CO2')
        co2.rotate(90, 'x')
        add_adsorbate(new_atoms, co2, h, position=position[0:2], mol_index=0)
    elif ads_name == 'H':
        H = Atoms('H')
        add_adsorbate(new_atoms, H, h, position=position[0:2])
    elif ads_name == 'H2O':
        h2o = molecule('H2O')
        h2o.rotate(180,'y')
        add_adsorbate(new_atoms, h2o, h, position=position[0:2])
    elif ads_name == 'OCCO':
        occo = read('/home/cat/vijays/tools/adsorbates/occo.xyz')
        add_adsorbate(new_atoms, occo, h, position=position[0:2])

    return new_atoms

def pickle2cube(homedir, atoms):
    import pickle
    import numpy as np
    from ase.io import read, write
    cd = pickle.load(open(homedir + '/density.pickle','rb'), encoding='latin1')

    nx,ny,nz = np.shape(cd)
    #cut away periodic image planes to correct QE output
    u=nx-1
    v=ny-1
    w=nz-1
    cd2 = np.empty((u,v,w), np.float)
    for i in range(u):
        for j in range(v):
            cd2[i][j][:] = cd[i][j][:w]
    write(homedir + '/density.cube',atoms,data=cd2)

    #edit density.cube grid size if odd number of grid points to correct for old versions of ASE
    bohr = 0.52917721092
    cell = atoms.get_cell()
    dx = cell[0][0]/(u*bohr)
    dy = cell[1][1]/(v*bohr)
    dz = cell[2][2]/(w*bohr)

    f = open(homedir + '/density.cube','r')
    lines = f.readlines()
    f.close()

    line3 = "%5.0f    %.6f    %.6f    %.6f\n"%(u,dx,0,0)
    line4 = "%5.0f    %.6f    %.6f    %.6f\n"%(v,0,dy,0)
    line5 = "%5.0f    %.6f    %.6f    %.6f\n"%(w,0,0,dz)

    lines[3] = line3
    lines[4] = line4
    lines[5] = line5

    f = open(homedir + '/density.cube','w')
    f.writelines(lines)
    f.close()

    return True

def find_max_empty_space(atoms, edir = 3):
    '''
    Assuming periodic boundary conditions, finds the largest
    continuous segment of free, unoccupied space and returns
    its midpoint in scaled coordinates (0 to 1) in the edir
    direction (default z)
    '''
    import numpy as np
    position_array = atoms.get_scaled_positions()[...,edir-1]##because index starts at 0
    position_array.sort()
    differences = np.diff(position_array)
    differences = np.append(differences, position_array[0] + 1-position_array[-1])##pbc
    max_diff_index = np.argmax(differences)
    if max_diff_index == len(position_array) -1:
        return (position_array[0]+1+position_array[-1])/2 %1
    else:
        return (position_array[max_diff_index] + position_array[max_diff_index+1])/2.


class AutoVivification(dict):
   def __getitem__(self, item):
       try:
           return dict.__getitem__(self, item)
       except KeyError:
           value = self[item] = type(self)()
           return value

def create_output_directory(output='output'):
    from pathlib import Path
    Path(output).mkdir(parents=True, exist_ok=True)
