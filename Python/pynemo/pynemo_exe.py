'''
Entry for the project

@author: Mr. Srikanth Nagella
'''

import sys, getopt
import profile
import logging

# Logging set to info
logging.basicConfig(level=logging.INFO)

def main():
    """ Main function which checks the command line parameters and
        passes them to the profile module for processing """

    setup_file = ''
    mask_gui = False
    try:
        opts, dummy_args = getopt.getopt(sys.argv[1:], "hgs:", ["mask_gui", "setup="])
    except getopt.GetoptError:
        print "usage: pynemo -s <namelist.bdy> "
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h":
            print "usage: pynemo [-g] -s <namelist.bdy> "
            sys.exit()
        elif opt in ("-s", "--setup"):
            setup_file = arg
        elif opt in("-g", "--mask_gui"):
            mask_gui = True

    if setup_file == "":
        print "usage: pynemo [-g] -s <namelist.bdy> "
        sys.exit(2)

    #Logger
    #logger = logging.getLogger(__name__)

    profile.go(setup_file, mask_gui)

if __name__ == "__main__":
    main()
