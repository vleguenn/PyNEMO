########################################################
## Parses a file to find nemo boundary settings to use #
#######                                          #######
## Written by John Farey, August 2012                 ##
########################################################
from __builtin__ import str

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
        self.filename = setfile
        if not setfile: # debug
            self.filename = '../data/namelist.bdy'
        try:
            namelist = open(self.filename, 'r')
        except:
            logger.error("Cannot open the file:"+setfile)
            raise
        data = namelist.readlines()
        
        # Dictionary of all the vars in the file and a seperate settings for boolean values        
        self.settings, self.bool_settings = self._assign(self._trim(data))
        
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
    def _get_val(self, vars, bool_vars, x):
        t = x[0][0:2].lower()
        n = x[0][3:].lower().strip() # 3 -> 0 to keep type info
        v = x[1].strip()

       
        if t == 'ln':
            if v.find('true') is not -1:
                if vars.has_key(n) != True:
                    vars[n] = True
                bool_vars[n] = True
            elif v.find('false') is not -1:
                if vars.has_key(n) != True:
                    vars[n] = False
                bool_vars[n] = False
            else:
                raise ValueError('Cannot assign %s to %s, must be boolean' 
                                                                    %(v,n))
        elif t == 'rn' or t == 'nn':
            if v.find('.') > -1 or v.find('e') > -1:
                try:
                    vars[n] = float(v)
                except ValueError:
                    print 'Cannot assign %s to %s, must be float'
                    raise
            else:
                try:
                    vars[n] = int(v)
                except ValueError:
                    print 'Cannot assign %s to %s, must be integer'
                    raise
        elif t == 'cn' or t == 'sn':              
                vars[n] = v.strip("'")                
        elif t == 'cl':
                name = x[0].split('(')
                index = name[1].split(')')
                if name[0].strip() not in vars.keys():
                    vars[name[0].strip()] = {}
                vars[name[0].strip()][index[0].strip()]=v.strip()
        else:
            raise ValueError('%s data type is ambiguous' %x)


    def _get_var_name_value(self,line):
        x = line.split("=",2)
        t = x[0][0:2].lower()
        name = x[0][3:].lower().strip() # 3 -> 0 to keep type info    
        v = x[1].strip()
        index = '-1'
        value = ''
        if t == 'ln':
            if v.find('true') is not -1:
                value = True
            elif v.find('false') is not -1:
                value = False
            else:
                raise ValueError('Cannot assign %s to %s, must be boolean' 
                                                                    %(v,name))
        elif t == 'rn' or t == 'nn':
            if v.find('.') > -1 or v.find('e') > -1:
                try:
                    value = float(v)
                except ValueError:
                    print 'Cannot assign %s to %s, must be float'
                    raise
            else:
                try:
                    value = int(v)
                except ValueError:
                    print 'Cannot assign %s to %s, must be integer'
                    raise
        elif t == 'cn' or t == 'sn':              
                value = v.strip("'")                
        elif t == 'cl':
                name = x[0].split('(')
                index = name[1].split(')')
                name = name[0].strip()
                index = index[0].strip()
                value = v                
        else:
            raise ValueError('%s data type is ambiguous' %x)
        return name, index, value

                
    def _replace_var_value(self,original_line,name,index,value, new_value):
        if type(value).__name__ == 'bool':
            value = str(value).lower()
            new_value = str(new_value).lower()
        return original_line.replace(value,new_value,1)
    
    def _get_var_name(self,line):
        x = line.split("=",2)
        t = x[0][0:2].lower()
        name = x[0][3:].lower().strip() # 3 -> 0 to keep type info        
        if t in ['ln', 'rn','nn','cn','sn']:
            return name, -1
        elif t == 'cl':
            name = x[0].split('(')
            index = name[1].split(')')
            return name, index
    
    # Returns tidy dictionary of var names and values
    def _assign(self, data):
        vars = {}
        bool_vars = {}
        for line in data:
            x = line.split('=', 2)
            self._get_val(vars, bool_vars, x)
#            n,v = self._get_val(x)
#            vars[n] = v

        return vars, bool_vars

    def strip_comments(self,line):
        line = line.rsplit('!')[0].strip()
        return line
    
        
    def write(self):
        '''
        This method write backs the variable data back into the file
        '''
        try:
            namelist = open(self.filename, 'r')
        except:
            self.logger.error("Cannot open the file:"+self.setfile+" to write ")
            raise
        data = namelist.readlines()
        
        for name in self.settings:
            values = self.settings[name]
            if type(values).__name__ != 'dict':
                values =  {'-1': values}
            for index, value in values.iteritems():
                count = -1
                var_found = False
                for line in data:
                    count = count + 1
                    #find the variable
                    line_without_comments = self.strip_comments(line)
                    if line_without_comments == '':
                        continue
                    print line_without_comments
                    data_name, data_index, data_value = self._get_var_name_value(line_without_comments)
                    
                    if data_name == name:
                        #found the variable line
                        if data_index == index:
                            if type(data_value).__name__ == 'bool' and type(value).__name__ != 'bool':
                                  data[count]=self._replace_var_value(line, name, index, data_value, self.bool_settings[name])
                                  continue
                            var_found = True
                            if data_value == value:
                                break
                            else:
                                data[count]=self._replace_var_value(line, name, index, data_value, value)
        namelist.close()
        namelist = open(self.filename,'w')
        namelist.truncate()
        namelist.writelines(data)
        namelist.close()                      
                #variable not found need to append to the file
#                if not var_found:
#                    data.append()
                                
                    #write back the new value
                

