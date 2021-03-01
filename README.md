ElectroParser
============

This package stores important quantities, like atoms objects into ASE databases. In addition, quantities such as surface dipole moments, workfunctions will also be stored. 

Note: These scripts were / are written during my PhD and are in continuous developement, feel free to clone, develop and make your own! ASE Databases make parsing calculations for plotting exceptionally easy, more details can be found here: https://wiki.fysik.dtu.dk/ase/ase/db/db.html

## Details

Three codes will be parsed 
* VASP
* QuantumEspresso
* GPAW

Quantities stored are:
1. Atoms: If you have 'OUTCAR', 'POSCAR', 'CONTCAR', 'vasprun.xml' it is a VASP calculation, if 'qn.traj', 'bfgs.traj', 'spe.traj', 'aiida.out', 'log' it is a QuantumEspresso calculation; if 'gpaw.traj', 'out.txt' it is a GPAW calculation
2. Workfunctions (if available)
3. Electric Field (if available) in $V/\AA$
4. Dipole moment in $e\AA$
5. Charge in $e$
6. Input file: Store the input file as a list of text 
7. paw: Plane wave cutoff 
8. Functional: Reads dispersion correction
9. ldau_data: Data for Hubbard-U calculations
10. Spin polarisation
11. d-band centre (deprecated)
12. Width of d-band (deprecated)
13. Hilbert transform (deprecated)
14. pdos: If there is a file called pdos.json, it is stored, helpful to store pdos important for your project


##  Usage

### Parser

The simplest way to put everything into an ASE database is to do the following
```
python -m parser --dbname <database-name> --levels <levels> --consider 'folder/name/here/'
```
The three options are:
1. dbname: Database name that you would like to store data in
2. levels: how many levels deep would you like the parser look into
3. consider (optional): Specific folder that you want to parse, the parser won't look elsewhere

### Restarts

Check if a *VASP* calculation has completed, if not tries to restart the calculation 

```
python -m check_restart <levels>
```

### Forces

Run the force method to get beta for electrochemical calculations from a single calculation
(Manuscript currently being written)

Call the ForceExtrapolation class
