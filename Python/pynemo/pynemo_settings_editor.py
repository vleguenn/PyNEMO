'''
Created on 7 Jan 2015

@author: Mr. Srikanth Nagella
'''

from PyQt4 import QtGui, QtCore
from gui.nemo_bdy_namelist_edit import NameListEditor
import nemo_bdy_setup

import sys
def main():
    app = QtGui.QApplication(sys.argv)
    fname = QtGui.QFileDialog.getOpenFileName(None, 'Open file')    
    setup = nemo_bdy_setup.Setup(fname)#'../../data/namelisttest.bdy')
    ex = NameListEditor(setup)
    ex.btnCancel.clicked.connect(lambda :sys.exit(0))
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    main()