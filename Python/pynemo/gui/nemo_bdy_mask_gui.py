'''
Created on 12 Jan 2015

@author: Mr. Srikanth Nagella
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
from .selection_editor import PolygonEditor, BoxEditor
import os.path
from PyQt4.QtCore import pyqtSignal, pyqtSlot, Qt
from nemo_bdy_mask import Mask
import logging
from PyQt4.QtGui import QSizePolicy
from matplotlib.colors import Normalize

mask_alpha = 0.3

from PyQt4 import QtGui
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.path import Path
from matplotlib.transforms import Bbox
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# pylint: disable=E1002
class MatplotlibWidget(QtGui.QWidget):
    """This class is a QWidget for pyNEMO mask plot"""
    min_depth = 200.0
    shelfbreak_dist = 200.0
    mask_type = 0
    def __init__(self, parent=None, mask=None, min_depth = 200.0, shelfbreak_dist = 200.0,*args, **kwargs):
        """ Initialises the mask, matplot and the navigation toolbar """
        super(MatplotlibWidget, self).__init__(parent)
        #QtGui.QWidget.__init__(self, parent)
        self.figure = Figure(*args, **kwargs)
        self.canvas = FigureCanvas(self.figure)
        self.mask = mask
        self.min_depth = min_depth
        self.shelfbreak_dist = shelfbreak_dist
        if self.mask is not None:
            self.mask.min_depth = min_depth
            self.mask.shelfbreak_dist = shelfbreak_dist
        self.toolbar = NemoNavigationToolbar(self.canvas, self)
        self.toolbar.locLabel.setMinimumWidth(100)
        self.toolbar.locLabel.setMaximumWidth(170)
        self.toolbar.locLabel.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
        self.toolbar.locLabel.setAlignment(Qt.AlignLeft|Qt.AlignTop)
        self.toolbar.drawing_tool.connect(self.drawing_tool_callback)
        self.axes = self.figure.add_subplot(111)
        self.cbar = None
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self._drawing_tool = None
        self._drawing_tool_name = None
        self.create_basemap()

    @pyqtSlot(str)
    def drawing_tool_callback(self, toolname):
        """ callback for the drawing tool when the signal of change of drawing tool is
        received"""
        if self._drawing_tool_name != None and toolname == "": #if tool is disabled
            self._drawing_tool.disable()
            self._drawing_tool_name = None
            self._drawing_tool = None
            self.canvas.draw()
        else:
            self._drawing_tool_name = toolname
            if self._drawing_tool_name == "freehand": #if freehand tool is enabled
                self._drawing_tool = PolygonEditor(self.axes, self.canvas)
                self.canvas.draw()
            elif self._drawing_tool_name == "rectangle": #if rectange tool is enabled
                self._drawing_tool = BoxEditor(self.axes, self.canvas)
                self._drawing_tool.enable()
                self.canvas.draw()

    def create_basemap(self):
        """ Draws the basemap and contour with mask information"""
        if self.mask == None:
            return

        x = np.arange(0, self.mask.lon.shape[0])
        y = np.arange(0, self.mask.lon.shape[1])
        x_vals, y_vals = np.meshgrid(y, x)
        Z = self.mask.bathy_data[...].astype(np.float64)
        #Z[Z==0] = np.nan
        Z = np.ma.masked_where(Z==0, Z)
        cmap = plt.get_cmap('GnBu')
        cmap.set_bad('0.0')
        cmap.set_under('black',1.0)
        cmap.set_over('black',1.0)
        transcmap = plt.get_cmap('autumn')
        transcmap.set_bad(alpha=0.5)
        masklayer = np.ma.masked_where(self.mask.data==-1,self.mask.data)
        extent = (0, self.mask.lon.shape[1],0, self.mask.lon.shape[0])
        #cax = self.axes.pcolormesh(x_vals, y_vals, Z, cmap=cmap)#, extend='min')#cmap=plt.get_cmap('GnBu'))#cmap=cm.s3pcpn)
        cax = self.axes.imshow(Z, cmap = cmap, origin="lower",extent=extent,aspect='auto')
        #self.axes.contourf(x_vals, y_vals, masklayer, [-2, -1, 0, 1, 2], cmap=transcmap,\
        #                   alpha=mask_alpha)
        self.axes.imshow(masklayer, cmap=transcmap,alpha=0.3,origin="lower",extent=extent,aspect='auto')

        zmin = np.amin(Z)
        zmax = np.amax(Z)
        if self.cbar is not None:
            self.cbar.remove()
        self.cbar = self.figure.colorbar(cax,ticks=np.linspace(zmin,zmax,10),orientation='horizontal')
        self.cbar.set_label("Bathymetry (units=%s)"%self.mask.data_units)
        self.canvas.draw()

    
    def reset_mask(self):
        if self.mask == None:
            return             
        self.mask.reset_mask()
        self.axes.clear()
        self.create_basemap()
        
    def add_mask(self):
        """ adds the selected region in the drawing tool to the mask """
        if self._drawing_tool_name != "" and self.mask != None:
            if self._drawing_tool.polygon != None:
                x = np.arange(0, self.mask.lon.shape[0])
                y = np.arange(0, self.mask.lon.shape[1])
                x_vals, y_vals = np.meshgrid(y, x)
                grid = zip(x_vals.ravel(), y_vals.ravel())
                
                self._drawing_tool.polygon.set_linewidth(1.0)
                p_path = Path(self._drawing_tool.polygon.xy)
                index = p_path.contains_points(grid)
                index = index.reshape(self.mask.lon.shape)
                xmin, ymin = np.min(self._drawing_tool.polygon.xy, axis=0)
                xmax, ymax = np.max(self._drawing_tool.polygon.xy, axis=0)
                self.mask.add_mask(index,[xmin,xmax,ymin,ymax])
                self._drawing_tool.reset()
                self.axes.clear()
                self.create_basemap()

    def remove_mask(self):
        """ removes the selected region in the drawing tool from the mask """
        if self._drawing_tool_name != "" and self.mask != None:
            if self._drawing_tool.polygon != None:
                x = np.arange(0, self.mask.lon.shape[0])
                y = np.arange(0, self.mask.lon.shape[1])
                x_vals, y_vals = np.meshgrid(y, x)
                grid = zip(x_vals.ravel(), y_vals.ravel()) #check for the index

                self._drawing_tool.polygon.set_linewidth(1.0)
                p_path = Path(self._drawing_tool.polygon.xy)
                index = p_path.contains_points(grid)
                index = index.reshape(self.mask.lon.shape)
                xmin, ymin = np.min(self._drawing_tool.polygon.xy, axis=0)
                xmax, ymax = np.max(self._drawing_tool.polygon.xy, axis=0)                
                self.mask.remove_mask(index,[xmin,xmax,ymin,ymax])
                self._drawing_tool.reset()
                self.axes.clear()
                self.create_basemap()

    def apply_border_mask(self):
        """ This applies an mask of given number of pixels at the border of the mask"""
        pixels, ok_btn_pressed = QtGui.QInputDialog.getText(self, 'Mask: Border Input',
                                                            'Enter number of pixel of border \
                                                             to be added to mask:')
        if ok_btn_pressed:
            self.mask.apply_border_mask(int(pixels))
            self.axes.clear()
            self.create_basemap()

    def set_mask_type(self,type):
        """ Sets the mask type """
        self.mask_type = type
        self.mask.mask_type = type
        
    @pyqtSlot(str, str)
    def set_bathymetry_file(self, bathymetry_filename, mask_file):
        """ Set the bathymetry file """
        try:
            self.mask = Mask(bathymetry_filename, mask_file, self.min_depth, self.shelfbreak_dist)
            self.mask.mask_type = self.mask_type
            self.create_basemap()
        except RuntimeError:
            pass # couldn't set the new file name
        
    @pyqtSlot(str)
    def save_mask_file(self, mask_file):
        """ Save the mask data to mask_file """
        if self.mask is not None:
            self.mask.save_mask(mask_file)
            
    @pyqtSlot(float, float)
    def set_mask_settings(self, min_depth, shelfbreak_dist):
        """ Mask settings update """
        self.min_depth = min_depth
        self.shelfbreak_dist = shelfbreak_dist
        self.mask.min_depth = min_depth
        self.mask.shelfbreak_dist = shelfbreak_dist
        
class NemoNavigationToolbar(NavigationToolbar):
    """ This is custom toolbar for the nemo which includes additional buttons
    for drawing tool and (add,remove) for mask in addtion to default NavigationToolbar
    provided by matplotlib """

    drawing_tool = pyqtSignal(str) #signal for the drawing tool changed
    def __init__(self, canvas, parent):
        """ Initialises the toolbar """
        self.toolitems = (('Home', 'Reset original view', 'home', 'home'),\
                          ('Back', 'Back to  previous view', 'back', 'back'),\
                          ('Forward', 'Forward to next view', 'forward', 'forward'),\
                          (None, None, None, None),\
                          ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),\
                          ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),\
                          ('Reset', 'Reset the mask', 'reset','reset'),\
                          (None, None, None, None),\
                          ('Freehand', 'Freehand drawing', 'freehand', 'freehand'),\
                          ('Rectangle', 'Rectangle drawing', 'rectangle', 'rectangle'),\
                          ('Border', 'Border selection', 'border', 'border'),\
                          ('plus', 'Add mask', 'add_mask', 'add_mask'),\
                          ('minus', 'Remove mask', 'remove_mask', 'remove_mask'),\
                          (None, None, None, None),\
                          ('Normal','Normal Mask','normal_mask','normal_mask'),\
                          ('MaxDepth', 'Max Depth Mask', 'max_depth_mask', 'max_depth_mask'),\
                          ('ShelfBreak','Shelf Break Mask','shelf_break_mask','shelf_break_mask'),\
                          (None, None, None, None)\
                          )
        NavigationToolbar.__init__(self, canvas, parent)
        self._actions['reset'].setIcon(set_icon('reset.png'))
        self._actions['freehand'].setCheckable(True)
        self._actions['freehand'].setIcon(set_icon('freehand.png'))
        self._actions['rectangle'].setCheckable(True)
        self._actions['rectangle'].setIcon(set_icon('rectangle.png'))
        self._actions['border'].setIcon(set_icon('border.png'))
        self._actions['add_mask'].setIcon(set_icon('plus.png'))
        self._actions['remove_mask'].setIcon(set_icon('minus.png'))
        self._actions['normal_mask'].setIcon((set_icon('all_mask.png')))
        self._actions['normal_mask'].setCheckable(True)
        self._actions['max_depth_mask'].setIcon((set_icon('max_depth.png')))
        self._actions['max_depth_mask'].setCheckable(True)
        self._actions['shelf_break_mask'].setIcon((set_icon('shelf_break.png')))
        self._actions['shelf_break_mask'].setCheckable(True)
        self.update_height_mask(0)
        
    def reset(self, *dummy):
        """ Callback for reset button clicked"""
        self.parent.reset_mask()

    def freehand(self, *dummy):
        """ callback for freehand button clicked """
        if self._actions['freehand'].isChecked() == True:
            if self._active == "PAN":
                self.pan()
            elif self._active == "ZOOM":
                self.zoom()
            elif self._actions['rectangle'].isChecked() == True:
                self._actions['rectangle'].setChecked(False)
                self.drawing_tool.emit("") # clear the rectangle selector
            self._active = None
            self.drawing_tool.emit('freehand')
            self._update_buttons_checked()
        else:
            self.drawing_tool.emit("")

    def rectangle(self, *dummy):
        """ callback for rectangel button clicked """
        if self._actions['rectangle'].isChecked() == True:
            if self._active == "PAN":
                self.pan()
            elif self._active == "ZOOM":
                self.zoom()
            elif self._actions['freehand'].isChecked() == True:
                self._actions['freehand'].setChecked(False)
                self.drawing_tool.emit("") # clear the freehand selector
            self._active = None
            self.drawing_tool.emit('rectangle')
            self._update_buttons_checked()
        else:
            self.drawing_tool.emit("")

    def border(self, *dummy):
        """ callback for border button clicked """
        self.parent.apply_border_mask()

    def add_mask(self, *dummy):
        """ callback for add mask button clicked """
        self.parent.add_mask()

    def remove_mask(self, *dummy):
        """ callback for remove mask button clicked """
        self.parent.remove_mask()

    def get_active_button(self):
        """ returns the current active button between freehand and rectangle"""
        if self._actions['rectangle'].isChecked() == True:
            return 'rectangle'
        elif self._actions['freehand'].isChecked() == True:
            return 'freehand'
        return None
    
    def normal_mask(self, *dummy):
        """ enable the normal mask button """
        self.update_height_mask(0)
    
    def max_depth_mask(self, *dummy):
        """ enables the minimum height mask """
        self.update_height_mask(1)
    
    def shelf_break_mask(self, *dummy):
        """ enables the shelf break mask button """
        self.update_height_mask(2)
    
    def update_height_mask(self, btn_id):
        """ update the height mask buttons in the interface """
        self._actions['normal_mask'].setChecked(False)
        self._actions['max_depth_mask'].setChecked(False)
        self._actions['shelf_break_mask'].setChecked(False)
        try:
            self.parent.set_mask_type(btn_id)
        except AttributeError:
            pass
        if btn_id == 0:
            self._actions['normal_mask'].setChecked(True)
        elif btn_id == 1:
            self._actions['max_depth_mask'].setChecked(True)
        elif btn_id == 2:
            self._actions['shelf_break_mask'].setChecked(True)

def set_icon(name):
    """ Creates an icon based on the file found in the module directory with input name"""
    return QtGui.QIcon(os.path.join(os.path.dirname(__file__), name))

