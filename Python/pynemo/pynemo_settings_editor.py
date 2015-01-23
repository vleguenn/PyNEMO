'''
Created on 7 Jan 2015

@author: Mr. Srikanth Nagella
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
from PyQt4 import QtGui

from .gui.nemo_bdy_input_window import InputWindow
import nemo_bdy_setup

import sys
def main():
    """ Main method which starts a Qt application and gives user
    an option to pick a namelist.bdy file to edit. Once user selects it
    it will open a dialog box where users can edit the parameters"""

    app = QtGui.QApplication(sys.argv)
    fname = QtGui.QFileDialog.getOpenFileName(None, 'Open file')
    setup = nemo_bdy_setup.Setup(fname)#'../../data/namelisttest.bdy')
    ex = InputWindow(setup)
    ex.nl_editor.btn_cancel.clicked.connect(lambda: sys.exit(0))
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
