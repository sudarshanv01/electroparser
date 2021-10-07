import sys, os
from electroparser.cli.parser_functions import data_to_store
from ase.db import connect
import click

@click.command()
@click.option('--dbname', default='test.db')
@click.option('--levels', default=1)
@click.option('--consider', default='')
def main(dbname, levels, consider):
    """Main function to store data in ASE database

    :param dbname: Name of database
    :type dbname: str
    :param levels: Maximum number of levels to look at
    :type levels: int
    :param consider: Each folder must contain this string
    :type consider: str
    """
    db = connect(dbname)
    all_data = data_to_store(levels, db, consider)



if __name__ == '__main__':
    main()



    
