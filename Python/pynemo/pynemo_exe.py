'''
Entry for the project

@author: Mr. Srikanth Nagella
'''

import sys, getopt
import pynemo.profile
import logging

# Logging set to info
logging.basicConfig(level=logging.INFO)

def main():
    """ Main function which checks the command line parameters and
        passes them to the profile module for processing """

    setup_file = ''
    try:
        opts = getopt.getopt(sys.argv[1:], "hs:", ["setup="])
    except getopt.GetoptError:
        print "usage: pynemo -s <namelist.bdy> "
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h":
            print "usage: pynemo -s <namelist.bdy> "
            sys.exit()
        elif opt in ("-s", "--setup"):
            setup_file = arg

    if setup_file == "":
        print "usage: pynemo -s <namelist.bdy> "
        sys.exit(2)

    #Logger
    #logger = logging.getLogger(__name__)

    pynemo.profile.go(setup_file)

if __name__ == "__main__":
    main()
