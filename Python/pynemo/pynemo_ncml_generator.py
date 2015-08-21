'''
Created on 2 Jul 2015

The main application object for hosting the pynemo ncml editor.
Used for development purposes to display the ncml editor dialog.

@author: Shirley Crompton, UK Science and Technology Facilities Council
'''
import sys
from PyQt4.QtGui import *
from gui import nemo_ncml_generator as ncml_generator
import logging
# Logging set to info
logging.basicConfig(level=logging.INFO)

def main():
    """ Command line execution method which check the input arguments and passes on to
    method to open the ncml generator window"""
    app = QApplication(sys.argv)
    ex = ncml_generator.Ncml_generator(None)
    sys.exit(app.exec_())

if __name__ == '__main__':
    main() 