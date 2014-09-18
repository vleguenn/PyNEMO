##
## written by John Farey, started 23 August 2012
##

'''
Reads a text file and deciphers variables.
It is a lot tidier just keeping all the variables from the file in a dictionary
rather than individually defining them

'''



#File to read
namelist = open('F:/NEMO_bdy_tools/first/namelist.bdy', 'r')

data = namelist.readlines()

namelist.close()

# Dictionary of all the vars in the file
#settings = {}



# Trims comments and whitespace
def _trim(data):
    for l in range(len(data)):
        x = data.pop(0)
        print 'x = ', x
        x = x.rsplit('!')[0].strip()
        if x is not '':
            data.append(x)

    return data


# x = string, checks type, returns name and val. Raises error for ambiguous
def _get_val(x):
    t = x[0][0:2].lower()
    n = x[0].lower()
    v = x[2]

    if t == 'ln':
        if v.find('true') is not -1:
            return n, True
        elif v.find('false') is not -1:
            return n, False
        else:
            raise ValueError('Cannot assign %s to %s, must be boolean', v,n)
    
    elif t is 'rn' or t is 'nn':
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

    elif t == 'cn' or 'sn':
        return n,v
    
    else:
        raise ValueError('"%s" is invalid assignment- data type is ambiguous', x)


# Returns tidy dictionary of var names and values
def _assign(data):
    vars = {}
    for l in data:
        x = l.split()
        n,v = _get_val(x)
        vars[n] = v

    return vars










settings = _assign(_trim(data))

