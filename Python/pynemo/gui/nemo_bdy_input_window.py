'''
Created on 21 Jan 2015

@author: Mr. Srikanth Nagella
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
# pylint: disable=E1002
from PyQt4 import QtGui
from .nemo_bdy_namelist_edit import NameListEditor
from .nemo_bdy_mask_gui import MatplotlibWidget

class InputWindow(QtGui.QDialog):
    '''
    Input Window for editing pyNEMO settings
    '''

    def __init__(self, setup):
        '''
        Initialises the UI components
        '''
        super(InputWindow, self).__init__()
        #initialise NameListEditor
        self.nl_editor = NameListEditor(setup)

        #initialise MatplotlibWidget
        self.mpl_widget = MatplotlibWidget()

        #connect namelistedit to matplotlibwidget
        self.nl_editor.mask_update.connect(self.mpl_widget.set_bathymetry_file)
        self.mpl_widget.set_bathymetry_file(setup.settings['bathy'])

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.nl_editor)
        hbox.addWidget(self.mpl_widget)

        self.setLayout(hbox)
        #set the Dialog title
        self.setWindowTitle("PyNEMO Settings Editor")
        #show the window
        self.show()
