#!/usr/bin/python

""" Check if a calculation has finished """
""" If it has not resubmit it """

from parser_functions import _get_homedir
from restart_functions import  _status_calc, _restart_calc
import sys, os


if __name__ == '__main__':

    level = int(sys.argv[1])
    if len(sys.argv) > 2:
        check_auto = sys.argv[2]
        if check_auto == 'T':
            autorestart = True
        else:
            autorestart = False
    else:
        autorestart = False
    user = os.environ.get('USER')

    homedirs = _get_homedir(level)

    for index, homedir in enumerate(homedirs):
        pwd = os.getcwd()
        try:
            status = _status_calc(homedir, user, pwd)
        except UnboundLocalError:
            continue

        print(homedir)

        if status == 'completed':
            print('Calculation has completed successfully')
        elif status == 'queue':
            print('Calculation in queue')
        elif status == 'failed':
            print('Calculation failed')
            if not autorestart:
                restart_flag = input('Should this calculation be restarted [Y/N]\n ')
                if restart_flag == 'Y':
                    _restart_calc(homedir)
                else:
                    print('Skipping')
            else:
                _restart_calc(homedir)
