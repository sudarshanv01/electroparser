ElectroParser
============

This package stores important quantities, like atoms objects into ASE databases. 

In addition, quantities such as surface dipole moments, workfunctions will also be stored. 

Three codes will be parsed 
* VASP
* QuantumEspresso
* GPAW

Note: These scripts were / are written during my PhD and are in continuous developement, feel free to clone, develop and make your own! ASE Databases make parsing calculations for plotting exceptionally easy, more details can be found here: https://wiki.fysik.dtu.dk/ase/ase/db/db.html

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
