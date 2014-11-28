'''
Entry for the project

@author: Mr. Srikanth Nagella
'''

import sys, getopt
import profile
import logging

# Logging set to info
logging.basicConfig(level=logging.INFO)

if __name__ == "__main__" :
    bathyfile = ''
    coordfile = ''
    setup_file = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:b:c:",["setup=","bathy=","coord="])
    except getopt.GetoptError:
        print "usage: pynemo -s <namelist.bdy> [-b <bathymetry file path> -c <output coord file path>]"
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == "-h":
            print "usage: pynemo -s <namelist.bdy> [-b <bathymetry file path> -c <output coord file path>]"
            sys.exit()
        elif opt in ("-s","--setup") :
            setup_file = arg
        elif opt in ("-b","--bathy") :
            bathyfile = arg
        elif opt in ("-c","--coord") :
            coordfile = arg
    if setup_file == "":
        print "usage: pynemo -s <namelist.bdy> [-b <bathymetry file path> -c <output coord file path>]"
        sys.exit(2)
        
    #Logger
    logger = logging.getLogger(__name__)
    if bathyfile == "":
        bathyfile = '../data/grid_C/NNA_R12_bathy_meter_bench.nc'
        logger.warn("No Bathymetry file provided. using ../data/grid_C/NNA_R12_bathy_meter_bench.nc")
    if coordfile == "":
        coordfile = '../data/toreador.nc'
        logger.warn("No output coordinate file provided. using ../data/toreador.nc")

    profile.go(setup_file,bathymeter_filepath=bathyfile,coord_output_filepath=coordfile)
