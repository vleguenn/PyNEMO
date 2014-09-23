########################################################
## Parses a file to find nemo boundary settings to use #
#######                                          #######
## Written by John Farey, August 2012                 ##
########################################################

'''
Invoke with a text file location, class init reads and deciphers variables.

attribute 'settings' is a dict holding all the vars.
'''

import os.path
import logging
import sys
class Setup:

    # setfile = settings file location
    def __init__(self, setfile):
        #Logging for class
        logger = logging.getLogger(__name__)
        
        if not setfile: # debug
            namelist = open('../data/namelist.bdy', 'r')
        else:
            try:
                namelist = open(setfile, 'r')
            except:
                logger.error("Cannot open the file:"+setfile)
                raise
        data = namelist.readlines()
        
        # Dictionary of all the vars in the file
        self.settings = self._assign(self._trim(data))

        namelist.close()

# # # # # # #
# Functions #
# # # # # # #

    # Trims comments and whitespace
    def _trim(self, data):
        for l in range(len(data)):
            x = data.pop(0)
            x = x.rsplit('!')[0].strip()
            if x is not '':
                data.append(x)

        return data

    # x = string, checks type, returns name and val. Raises error for ambiguous
    def _get_val(self, x):
        t = x[0][0:2].lower()
        n = x[0][3:].lower().strip() # 3 -> 0 to keep type info
        v = x[1].strip()

        if t == 'ln':
            if v.find('true') is not -1:
                return n, True
            elif v.find('false') is not -1:
                return n, False
            else:
                raise ValueError('Cannot assign %s to %s, must be boolean' 
                                                                    %(v,n))
        elif t == 'rn' or t == 'nn':
            if v.find('.') > -1 or v.find('e') > -1:
                try:
                    return n, float(v)
                except ValueError:
                    print 'Cannot assign %s to %s, must be float'
                    raise
            else:
                try:
                    return n, int(v)
                except ValueError:
                    print 'Cannot assign %s to %s, must be integer'
                    raise
        elif t == 'cn' or t == 'sn':
            return n,v.strip("'")
        else:
            raise ValueError('%s data type is ambiguous' %x)

    # Returns tidy dictionary of var names and values
    def _assign(self, data):
        vars = {}
        for line in data:
            x = line.split('=', 2)
            n,v = self._get_val(x)
            vars[n] = v

        return vars




