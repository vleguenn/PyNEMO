'''
Created on 2 Jul 2015

@author: Shirley Crompton, UK Science and Technology Facilities Council
'''
import logging
import os
from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtCore import pyqtSlot

class Ncml_tab(QtGui.QWidget):
    '''
    tab contents to define child aggregation
    '''
    def __init__(self, tabName):
        '''
        Initialises the UI components
        '''
        super(Ncml_tab, self).__init__() 
        self.logger = logging.getLogger(__name__)
        self.var= tabName  # tabName is used to determine the set of netcdf variables to process
        #print self.var
        #self.value_dict = {} # to capture the user inputs
    
        # no params yet, may be allow user to predefine an input ncml for edit????    
        self.initUI()
        
        
    def initUI(self):  
        QtGui.QToolTip.setFont(QtGui.QFont('SansSerif', 11))
        self.varStackedWidget = QtGui.QStackedWidget()                   
        #variable chooser combobox
        combo_vars = []
        if(self.var == unicode("Tracer").encode('utf-8')):
            combo_vars = [unicode('temperature').encode('utf-8'),unicode('salinity').encode('utf-8')] #votemper, vosaline
            self.votemper = ncml_variable(unicode('temperature').encode('utf-8'))
            self.vosaline = ncml_variable(unicode('salinity').encode('utf-8'))
            self.varStackedWidget.addWidget(self._addStackWidget())
            self.varStackedWidget.addWidget(self._addStackWidget())
            #debug
            print 'Tracer has ' + str(self.varStackedWidget.count())
        elif(self.var == unicode("Ice").encode('utf-8')):
            combo_vars = [unicode('ice thickness').encode('utf-8'),unicode('leads fraction').encode('utf-8'),unicode('snow thickness').encode('utf-8')] #'iicethic,ileadfra,isnowthi
            self.iicethic = ncml_variable(unicode('ice_thickness').encode('utf-8'))
            self.ileadfra = ncml_variable(unicode('leads_fraction').encode('utf-8'))
            self.isnowthi = ncml_variable(unicode('snow_thickness').encode('utf-8'))
            self.varStackedWidget.addWidget(self._addStackWidget())
            self.varStackedWidget.addWidget(self._addStackWidget())
            self.varStackedWidget.addWidget(self._addStackWidget())
            print 'Ice has ' + str(self.varStackedWidget.count())
        elif(self.var == unicode("Dynamics").encode('utf-8')):
            combo_vars = [unicode('zonal velocity').encode('utf-8'), unicode('meridian velocity').encode('utf-8'), unicode('sea surface height').encode('utf-8')] #vozocrtx, vomecrty, sossheig
            self.vozocrtx = ncml_variable(unicode('zonal_velocity').encode('utf-8'))
            self.vomecrty = ncml_variable(unicode('meridian_velocity').encode('utf-8'))
            self.sossheig = ncml_variable(unicode('sea_surface_height').encode('utf-8'))
            self.varStackedWidget.addWidget(self._addStackWidget())
            self.varStackedWidget.addWidget(self._addStackWidget())
            self.varStackedWidget.addWidget(self._addStackWidget())
            print 'Dynamics has ' + str(self.varStackedWidget.count())
        self.varStackedWidget.setCurrentIndex(0)  #we rely on the stacked tab index to be the same as the combo box 
        '''
        these two tabs are disabled
        elif(self.var == "Ecosysem"):
            vars = ['nitrate','silicate'] #nitrate, silicate
        elif(self.var == "Grid"):
        '''
        #combo box     
        self.var_combo = QtGui.QComboBox()
        self.var_combo.addItems(combo_vars)
        self.var_combo.setEditable(False)
        self.var_combo.setCurrentIndex(0)
        #the value if not saved is cached during the session, we can wait until the add button is pressed
        self.var_combo.currentIndexChanged.connect(lambda var_name = self.var : self.src_combo_changed(var_name))
        self.var_combo.currentIndexChanged.connect(self.setWidgetStack)
        #label
        var_label = QtGui.QLabel(unicode('Variable').encode('utf-8'))
        #set layout
        stacked_hBox = QtGui.QHBoxLayout()
        stacked_hBox.setMargin(5)
        stacked_hBox.setSpacing(50) # spacing between items
        stacked_hBox.setAlignment(QtCore.Qt.AlignLeft)
        stacked_hBox.addWidget(var_label)
        stacked_hBox.addWidget(self.var_combo)
        #
        vBoxLayout = QtGui.QVBoxLayout()        
        vBoxLayout.addLayout(stacked_hBox)
        vBoxLayout.addWidget(self.varStackedWidget)
        #
        grp_box = QtGui.QGroupBox(None)
        grp_box.setLayout(vBoxLayout)
                        
        '''
        :TODO Need to add the override time gui widgets
        '''
        
        
        '''
        button bar
        '''
        # reset button
        reset_btn = QtGui.QPushButton(unicode('Reset').encode('utf-8'))
        reset_btn.setToolTip(unicode('Reset fields to previously saved values').encode('utf-8'))
        add_btn = QtGui.QPushButton(unicode('Add').encode('utf-8')) 
        add_btn.setDefault(True)
        add_btn.setToolTip(unicode('Add the current definition to the NcML').encode('utf-8'))
        #connect up with events
        reset_btn.clicked.connect(self.reset_tab)
        add_btn.clicked.connect(self.add_tab)
       
        btn_hBox = QtGui.QHBoxLayout(None)
        btn_hBox.setMargin(5)
        btn_hBox.setSpacing(10)
        btn_hBox.setAlignment(QtCore.Qt.AlignCenter)
        btn_hBox.addWidget(reset_btn)
        btn_hBox.addWidget(add_btn)
        
        #build the contents         
        vbox = QtGui.QVBoxLayout(self)
        vbox.setSpacing(10)
        vbox.setContentsMargins(10, 10, 5, 5)
        vbox.addWidget(grp_box)
        vbox.addLayout(btn_hBox)
    '''
    create the stacked widget for each nemo variable
    '''    
    def _addStackWidget(self):
        self.varWidget = QtGui.QWidget()
        #self.varWidget.setObjectName(objName)
        varLayout = QtGui.QGridLayout()
        varLayout.setSpacing(20)
               
        #labels
        src_label = QtGui.QLabel(unicode('Source directory*').encode('utf-8'))  
        cbox_label = QtGui.QLabel(unicode('Includes subdirs').encode('utf-8'))
        regex_label = QtGui.QLabel(unicode('Regular expression').encode('utf-8'))
        old_name_label = QtGui.QLabel(unicode('Existing variable name*').encode('utf-8'))        
        #input textboxs
        self.varWidget.src_tedit = QtGui.QLineEdit()       # input widgets need to be attached to the stacked widget itself 
        self.varWidget.src_tedit.setToolTip(unicode('either remote OPeNDAP server or local file absolute path').encode('utf-8'))
        self.varWidget.src_tedit.returnPressed.connect(self.src_tedit_edited)
        
        
        self.varWidget.cbox = QtGui.QCheckBox()
        self.varWidget.cbox.setCheckable(True)
        self.varWidget.cbox.setChecked(False)
        self.varWidget.cbox.setToolTip(unicode('includes subdirectories').encode('utf-8'))
        self.varWidget.regex_tedit = QtGui.QLineEdit()
        self.varWidget.regex_tedit.setToolTip(unicode('see http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/ncml/AnnotatedSchema4.html#regexp').encode('utf-8'))                
        self.varWidget.old_name_tedit = QtGui.QLineEdit()
        self.varWidget.old_name_tedit.setToolTip(unicode('variable name in data file').encode('utf-8'))
        
        varLayout.addWidget(src_label, 1, 0, 1, 1)
        varLayout.addWidget(self.varWidget.src_tedit, 1, 1, 1, 3)
        varLayout.addWidget(cbox_label, 2, 0, 1, 1)
        varLayout.addWidget(self.varWidget.cbox, 2, 1, 1, 1)        
        varLayout.addWidget(regex_label, 2, 2, 1, 1)
        varLayout.addWidget(self.varWidget.regex_tedit, 2, 3, 1, 1)
        varLayout.addWidget(old_name_label, 3, 0, 1, 1)
        varLayout.addWidget(self.varWidget.old_name_tedit, 3, 1, 1, 3)
        
        self.varWidget.setLayout(varLayout)
        return self.varWidget
    '''
    synchronise stack widget display with combo box value changed callback
    '''
    @pyqtSlot()
    def setWidgetStack(self):
        self.varStackedWidget.setCurrentIndex(self.var_combo.currentIndex())
    '''
    variable combo box value changed callback
    '''
    @pyqtSlot()
    def src_combo_changed(self, var_name):  
        #not sure why the current text is prefixed by the index : eg 0temperature      
        print 'src_combo_value_changed to : ' + str(var_name) +  unicode(str(self.var_combo.currentText())).encode('utf_8')
        
        
    @pyqtSlot()
    def src_tedit_edited(self):
        src_tedit_input = self.varStackedWidget.currentWidget().src_tedit.text()
        print 'src_edit text edited : ', src_tedit_input
        #validate the input now
        if not str(src_tedit_input).startswith('http'): 
            if not os.path.isabs(src_tedit_input): #assumes local file
                QtGui.QMessageBox.critical(self, unicode('Something is wrong').encode('utf-8'), unicode('source directory must be an absolute path!').encode('utf-8'), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                self.varStackedWidget.currentWidget().src_tedit.clear()
            if not os.path.exists(src_tedit_input) :
                QtGui.QMessageBox.critical(self, unicode('Something is wrong').encode('utf-8'), unicode('source directory does not exist!').encode('utf-8'), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                self.varStackedWidget.currentWidget().src_tedit.clear()   
    '''
    reset button pushed callback.  The widgets are reset to default values
    '''
    @pyqtSlot()  
    def reset_tab(self):        
        #current screen value is not saved until the add button is pressed
        #reset only reset the screen values, not the cached values        
        if self.var_combo.currentText() == unicode("temperature").encode('utf-8'):
            print 'reset button is pushed, temperature ....'
            self.resetValues(self.votemper)
        elif self.var_combo.currentText() == unicode("salinity").encode('utf-8'):
            self.resetValues(self.vosaline)
        elif self.var_combo.currentText() == unicode("ice thickness").encode('utf-8'):
            self.resetValues(self.iicethic)
        elif self.var_combo.currentText() == unicode("leads fraction").encode('utf-8'):
            self.resetValues(self.ileadfra)
        elif self.var_combo.currentText() == unicode("snow thickness").encode('utf-8'):
            self.resetValues(self.isnowthi)
        elif self.var_combo.currentText() == unicode("zonal velocity").encode('utf-8'):
            self.resetValues(self.vozocrtx) 
        elif self.var_combo.currentText() == unicode("meridian velocity").encode('utf-8'):
            self.resetValues(self.vomecrty)
        elif self.var_combo.currentText() == unicode("sea surface height").encode('utf-8'):
            self.resetValues(self.sossheig)   
            
    '''
    reset the stacked widget values
    '''    
    def resetValues(self, currentValues = None):
        #print 'in resetValues ....'
        if currentValues is None:
            #self.var_combo.setCurrentIndex(0)    #we don't reset this, as this is the key
            self.varStackedWidget.currentWidget().src_tedit.clear()
            self.varStackedWidget.currentWidget().regex_tedit.clear()
            self.varStackedWidget.currentWidget().old_name_tedit.clear()
            self.varStackedWidget.currentWidget().cbox.setChecked(False)
        else:
            #print 'name : ' + currentValues.name + ', src: ' + currentValues.src + ', regex: ' +  currentValues.regex + ', old_name: ' + currentValues.old_name
            #self.var_combo.setCurrentIndex(0)
            self.varStackedWidget.currentWidget().src_tedit.setText(currentValues.src)
            self.varStackedWidget.currentWidget().regex_tedit.setText(currentValues.regex)
            self.varStackedWidget.currentWidget().old_name_tedit.setText(currentValues.old_name)
            self.varStackedWidget.currentWidget().cbox.setChecked(currentValues.subdirs)
        
    '''
    add button pushed call back
    '''
    @pyqtSlot()  
    def add_tab(self): 
        #first validate the src tab is not null
        if(self.varStackedWidget.currentWidget().src_tedit.text() is None or self.varStackedWidget.currentWidget().src_tedit.text() == '' or               
           self.varStackedWidget.currentWidget().old_name_tedit.text() is None or self.varStackedWidget.currentWidget().old_name_tedit.text() == ''):
                QtGui.QMessageBox.critical(self, unicode('Something is wrong').encode('utf-8'), unicode('source directory and existing variable name cannot be blank!').encode('utf-8'), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
        else:
            '''        
            if not str(target).startsWith(unicode('http').encode('utf-8')):
                if os.path.exists(os.path.normpath(target)) == False:
                    QMessageBox.critical(self, unicode('Something is wrong').encode('utf-8'), unicode('source directory does not exist!').encode('utf-8'), QMessageBox.Ok, QMessageBox.Ok)
                    return #breakout now 
            '''
            #print type(self.var) = str   
            if(self.var == unicode("Tracer").encode('utf-8')):                
                if (self.var_combo.currentText() == unicode("temperature").encode('utf-8')):            
                    if(self._sameValues(self.votemper)):
                        QtGui.QMessageBox.information(self, 'For information', 'No changes have been made!', QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else:   
                        self.votemper.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.votemper.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        self.votemper.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.votemper.regex = self.varStackedWidget.currentWidget().regex_tedit.text()
                        else:
                            self.votemper.regex = ''    #blank it over    
                else: # can only be salinity
                    if(self._sameValues(self.vosaline)):
                        QtGui.QMessageBox.information(self, 'For information', 'No changes have been made!', QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else:
                        self.vosaline.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.vosaline.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        self.vosaline.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.vosaline.regex = self.varStackedWidget.currentWidget().regex_tedit.text()                            
                        else:
                            self.vosaline.regex = ''
            elif(self.var == unicode('Ice').encode('utf-8')): #iicethic,ileadfra,isnowthi
                if (self.var_combo.currentText() == unicode("ice thickness").encode('utf-8')):            
                    if(self._sameValues(self.iicethic)):
                        QtGui.QMessageBox.information(self, 'For information', 'No changes have been made!', QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else:
                        self.iicethic.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.iicethic.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        self.iicethic.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.iicethic.regex = self.varStackedWidget.currentWidget().regex_tedit.text() 
                        else:
                            self.iicethic.regex = ''                       
                elif(self.var_combo.currentText() == unicode("leads fraction").encode('utf-8')): 
                    if(self._sameValues(self.ileadfra)):
                        QtGui.QMessageBox.information(self, unicode('For information').encode('utf-8'), unicode('No changes have been made!').encode('utf-8'), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else:
                        self.ileadfra.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.ileadfra.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        self.ileadfra.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.ileadfra.regex = self.varStackedWidget.currentWidget().regex_tedit.text()                            
                        else:
                            self.ileadfra.regex = ''
                else:
                    if(self._sameValues(self.isnowthi)): #snow thickness
                        QtGui.QMessageBox.information(self, unicode('For information').encode('utf-8'), unicode('No changes have been made!').encode('utf-8'), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else:
                        self.isnowthi.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.isnowthi.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        self.isnowthi.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.isnowthi.regex = self.varStackedWidget.currentWidget().regex_tedit.text()                #Dynamics ['zonal velocity', 'meridian velocity', 'sea surface height'] #vozocrtx, vomecrty, sossheig
                        else:
                            self.isnowthi.regex = ''
            elif(self.var == unicode("Dynamics").encode('utf-8')):
                if (self.var_combo.currentText() == unicode("zonal velocity").encode('utf-8')):            
                    if(self._sameValues(self.vozocrtx)):
                        QtGui.QMessageBox.information(self, unicode('For information').encode('utf-8'), unicode('No changes have been made!').encode('utf-8'), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else:
                        self.vozocrtx.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.vozocrtx.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        self.vozocrtx.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.vozocrtx.regex = self.varStackedWidget.currentWidget().regex_tedit.text()
                        else:
                            self.vozocrtx.regex = ''
                elif(self.var_combo.currentText() == unicode('meridian velocity').encode('utf-8')): 
                    if(self._sameValues(self.vomecrty)):
                        QtGui.QMessageBox.information(self, 'For information', 'No changes have been made!', QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else: 
                        self.vomecrty.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.vomecrty.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        self.vomecrty.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.vomecrty.regex = self.varStackedWidget.currentWidget().regex_tedit.text()
                        else:
                            self.vomecrty.regex = ''
                elif(self.var_combo.currentText() == unicode('sea surface height').encode('utf-8')):      
                    if(self._sameValues(self.sossheig)): #sea surface height
                        QtGui.QMessageBox.information(self, unicode('For information').encode('utf-8'), unicode('No changes have been made!').encode('utf-8'), QtGui.QMessageBox.Ok, QtGui.QMessageBox.Ok)
                    else:
                        self.sossheig.src = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
                        self.sossheig.subdirs = self.varStackedWidget.currentWidget().cbox.isChecked()
                        self.sossheig.old_name = self.varStackedWidget.currentWidget().old_name_tedit.text()
                        if(self.varStackedWidget.currentWidget().regex_tedit.text() is not None or self.varStackedWidget.currentWidget().regex_tedit.text() != ''):
                            self.sossheig.regex = self.varStackedWidget.currentWidget().regex_tedit.text()  
                        else:
                            self.sossheig.regex = ''
        
    '''
    convert the target folder into NcMl required format
    '''
    def _convertSrc(self, thepath):
        #print 'thepath before trimming ', thepath
        fpath = str(thepath.trimmed()) #24Aug15 trimmed whitespaces at both end of QString and then cast to string
        #print 'thepath after trimming and casting ', fpath
        #make sure thatit is an absolute path and prefixed with file:/ and uses / file separator
        if fpath.startswith('http:/'):
            target = fpath # do nothing 
        elif fpath.startswith('file:/'):
            temp = os.path.normpath(fpath[6:])
            print 'normal path : ', temp
            target =  unicode('file:/' + str(os.path.abspath(temp)).replace("\\", "/")).encode('utf-8') 
        else: #should be local file but not prefixed by file:/  we still check for absolute path
            print 'normal path : ', os.path.normpath(fpath)
            target = unicode('file:/' + str(os.path.abspath(fpath)).replace("\\", "/")).encode('utf-8')
        
        if not str(target).endswith('/'):
            target = target + '/'
            
        return target    
    
    '''   
    compare the gui cached values with the stored values
    '''   
    def _sameValues(self, ncml_var):
        print 'before state - variable: ' + ncml_var.name + ', src: ' + ncml_var.src + ', regex: ' +  ncml_var.regex + ', old_name: ' + ncml_var.old_name
        target = self._convertSrc(self.varStackedWidget.currentWidget().src_tedit.text())
        
        if(target == ncml_var.src and \
           self.varStackedWidget.currentWidget().old_name_tedit.text() is not None and self.varStackedWidget.currentWidget().old_name_tedit.text() == ncml_var.old_name and \
           self.varStackedWidget.currentWidget().regex_tedit.text() is not None and self.varStackedWidget.currentWidget().regex_tedit.text() == ncml_var.regex and \
           self.varStackedWidget.currentWidget().cbox.isChecked() == ncml_var.subdirs):
            return True
        else:
            return False
'''
convenient class to hold the user input for each variable
'''   
class ncml_variable(object):
    '''
    convenient class to hold the values for a ncml variable
    '''
    def __init__(self, varName):
        #print 'created ncml_variable object : ' + varName
        self.name = varName
        self.src = ''
        self.regex = ''
        self.old_name = ''
        self.subdirs = False
        
        
        
        
        