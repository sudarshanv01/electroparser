#!/usr/bin/python

import sys, os
from parser_functions import data_to_store
from ase.db import connect

if __name__ == '__main__':

    dbname = str(sys.argv[1])
    levels_max = int(sys.argv[2])

    db = connect(dbname)

    all_data = data_to_store(levels_max, db)

    
