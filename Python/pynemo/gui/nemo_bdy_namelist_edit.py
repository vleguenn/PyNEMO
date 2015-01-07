'''
Editor for namelist.bdy file

@author: Mr. Srikanth Nagella
'''

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import SIGNAL,SLOT 

import sys
import ast

class NameListEditor(QtGui.QDialog):
    '''
    This class creates a gui for the Namelist file options
    '''
    new_settings={}

    def __init__(self, setup):
        '''
        Constructor for setting up the gui using the settings
        '''
        super(NameListEditor, self).__init__()
        self.settings = setup.settings
        self.setup    = setup
        self.initUI()
        
    def initUI(self):
        '''
        Initialises the UI components of the GUI
        '''
        client = QtGui.QWidget(self)
        # Create the Layout to Grid
        grid = QtGui.QGridLayout()      
        
        # Loop through the settings and create widgets for each setting
        index = 0
        for setting in self.settings:
            # initialises setting Widget
            label = QtGui.QLabel(setting)
            if type(self.settings[setting]).__name__ in ['str','float','double','int','time','dict']:                                               
                text  = QtGui.QLineEdit(self)
                text.setText(str(self.settings[setting]))
                text.textChanged.connect(lambda value=setting,var_name=setting: self.labelChanged(value,var_name))
            elif type(self.settings[setting]).__name__ == 'bool':
                text = QtGui.QComboBox(self)
                text.insertItem(0,'True')
                text.insertItem(1,'False')                
                if self.settings[setting]:
                    text.setCurrentIndex(0)
                else:
                    text.setCurrentIndex(1)
                text.currentIndexChanged.connect(lambda value=setting, var_name=setting: self.comboIndexChanged(value,var_name))
            
            grid.addWidget(label, index, 0)
            grid.addWidget(text,  index, 1)
            index=index+1
            
        client.setLayout(grid)
        
        #scrollbars
        scrollArea = QtGui.QScrollArea(self)
        scrollArea.setWidget(client)
        
        #save cancel buttons
        btnWidget = QtGui.QWidget(self)
        hboxLayout = QtGui.QHBoxLayout(self)
        btnSave = QtGui.QPushButton('Save')
        btnSave.clicked.connect(self._btn_save_callback)
        self.btnCancel = QtGui.QPushButton('Cancel')
        self.btnCancel.clicked.connect(self._btn_cancel_callback)
        hboxLayout.addWidget(btnSave)
        hboxLayout.addWidget(self.btnCancel)
        btnWidget.setLayout(hboxLayout)
        
        boxLayout = QtGui.QVBoxLayout(self)
        boxLayout.addWidget(scrollArea)
        boxLayout.addWidget(btnWidget)
        
        
        self.setLayout(boxLayout) 
        #set the Dialog title
        self.setWindowTitle("PyNEMO Settings Editor")           
        #show the window           
        self.show()
        
    def labelChanged(self,value,name):
        self.new_settings[name] = unicode(value).encode('utf_8')
#        print self.new_settings
        
    def comboIndexChanged(self,value,name):
        if value == 0:
            self.new_settings[name] = True
        else:
            self.new_settings[name] = False
            
    def _btn_save_callback(self):
        #copy the the modified values to settings and call the setup save
        print self.settings
        for setting in self.new_settings:
            if (type(self.settings[setting]).__name__ == 'dict') & (type(self.new_settings[setting]).__name__ != 'dict') :
                self.new_settings[setting] = ast.literal_eval(self.new_settings[setting])
            self.settings[setting] = self.new_settings[setting]
        print self.settings
        self.setup.settings = self.settings
        self.setup.write()
    
    def _btn_cancel_callback(self):
        self.close()

