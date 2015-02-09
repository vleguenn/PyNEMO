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
from PyQt4.QtCore import pyqtSignal, pyqtSlot
import logging

mask_alpha = 0.3
class Mask(object):
    """This is a Mask holder. which reads from a netCDF bathemetry file and
    stores it in 'data' member variable"""

    data = None
    bathy_data = None
    mask_file = None
    def __init__(self, bathymetry_file=None, mask_file=None):
        """Initialises the Mask data"""
        self.set_mask_file(mask_file)
        self.set_bathymetry_file(bathymetry_file)
        self.logger = logging.getLogger(__name__)

    def set_mask_file(self, mask_file):
        """Reads the mask data from the mask file"""
        self.mask_file = mask_file
        #if mask file is not set then reset the data
        if self.mask_file == None:
            self.data = None
            return

        try:
            mask_nc = Dataset(str(self.mask_file), mode="r")
            data = mask_nc.variables['mask']
            self.data = data[:,:]
        except (IOError, RuntimeError):
            self.logger.error('Cannot open mask file')
            self.data = None

    def set_bathymetry_file(self, bathy_file):
        """ This reads the bathymetry file and sets the land to 0 and ocean to 1 """
        if bathy_file == None:
            return

        try:
            self.bathymetry_file = str(bathy_file)
            #open the bathymetry file
            self.bathy_nc = Dataset(self.bathymetry_file)
            self.lon = np.asarray(self.bathy_nc.variables['nav_lon'])
            self.lat = np.asarray(self.bathy_nc.variables['nav_lat'])
            if self.data is None:
                self.data = self.bathy_nc.variables['Bathymetry']
                self.data = np.asarray(self.data[:, :])
                self.data = np.around((self.data + .5).clip(0, 1))
            self.bathy_data = self.bathy_nc.variables['Bathymetry']
        except IOError:
            self.logger.error('Cannot open bathymetry file')

    def save_mask(self, mask_file):
        """Reads the mask data from the mask file"""
        if mask_file == None:
            mask_file = self.mask_file

        try:
            self.logger.info('Writing mask data to %s' % mask_file)
            mask_nc = Dataset(str(mask_file), mode="w")
            mask_nc.createDimension('y', size=len(self.bathy_nc.dimensions['y']))
            mask_nc.createDimension('x', size=len(self.bathy_nc.dimensions['x']))
            nav_lat = mask_nc.createVariable('nav_lat', 'f4', ('y','x',))
            nav_lon = mask_nc.createVariable('nav_lon', 'f4', ('y','x',))
            mask_var = mask_nc.createVariable('mask', 'f4', ('y','x',))            
            mask_var[...] = self.data
            nav_lat[...] = self.lat
            nav_lon[...] = self.lon
        except IOError:
            self.logger.info('Cannot open mask file for writing')

    def apply_border_mask(self, pixels):
        """ pixels is number of pixels in the border that need applying mask"""
        if self.data is not None and pixels < self.data.shape[0] and pixels < self.data.shape[1]:
            tmp = np.ones(self.data.shape, dtype=bool)
            tmp[pixels:-pixels, pixels:-pixels] = False
            self.data[tmp] = -1


    def apply_mediterrian_mask(self):
        """ This is mediterrian mask specific for the test bathymetry file """
        tmp = self.data[0:59, 280:350]
        tmp[tmp == 1] = -1
        self.data[0:59, 280:350] = tmp

from PyQt4 import QtGui
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.path import Path
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
# pylint: disable=E1002
class MatplotlibWidget(QtGui.QWidget):
    """This class is a QWidget for pyNEMO mask plot"""
    def __init__(self, parent=None, mask=None, *args, **kwargs):
        """ Initialises the mask, matplot and the navigation toolbar """
        super(MatplotlibWidget, self).__init__(parent)
        #QtGui.QWidget.__init__(self, parent)
        self.figure = Figure(*args, **kwargs)
        self.canvas = FigureCanvas(self.figure)
        self.mask = mask
        self.toolbar = NemoNavigationToolbar(self.canvas, self)
        self.toolbar.drawing_tool.connect(self.drawing_tool_callback)
        self.axes = self.figure.add_subplot(111)
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

#     def create_basemap(self):
#         """ Draws the basemap and contour with mask information"""
#         if self.mask == None:
#             return
#         #give some extra space for the mask to be shown on map
#         minlon = self.mask.lon.min()-1
#         maxlon = self.mask.lon.max()+1
#         minlat = self.mask.lat.min()-1
#         maxlat = self.mask.lat.max()+1
# #        basemap = Basemap(projection='cyl', llcrnrlat=minlat, urcrnrlat=maxlat,\
# #                    llcrnrlon=minlon, urcrnrlon=maxlon, resolution='c', ax=self.axes)
#
#         basemap = Basemap(projection='cea', llcrnrlat=minlat, urcrnrlat=maxlat,\
#                     llcrnrlon=minlon, urcrnrlon=maxlon, resolution='c', ax=self.axes)
#
# #        basemap = Basemap(projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,\
# #                    llcrnrlon=minlon, urcrnrlon=maxlon, resolution='c', ax=self.axes)
#         basemap.drawcoastlines()
#         # draw parallels and meridians.
#         basemap.drawparallels(np.arange(-90., 91., 30.))
#         basemap.drawmeridians(np.arange(-180., 181., 60.))
#         basemap.drawmapboundary(fill_color='aqua')
#         clevs = [-1, 0, 1]
#         bathy_clevs = np.arange(-1,5300,100)
#         x_vals, y_vals = basemap(self.mask.lon, self.mask.lat)
#         basemap.contourf(x_vals, y_vals, self.mask.bathy_data, bathy_clevs ,
#                            cmap=plt.get_cmap('GnBu'))#cmap=cm.s3pcpn)
#         basemap.contourf(x_vals, y_vals, self.mask.data, clevs,
#                            cmap=plt.get_cmap('autumn'), alpha=mask_alpha)
# #        basemap.contourf(x_vals, y_vals, )
# #        plt.title("Equidistant Cylindrical Projection")
#         self.canvas.draw()
    def create_basemap(self):
        """ Draws the basemap and contour with mask information"""
        if self.mask == None:
            return
        #give some extra space for the mask to be shown on map
        minlon = self.mask.lon.min()-1
        maxlon = self.mask.lon.max()+1
        minlat = self.mask.lat.min()-1
        maxlat = self.mask.lat.max()+1

#        basemap = Basemap(projection='cyl', llcrnrlat=minlat, urcrnrlat=maxlat,\
#                    llcrnrlon=minlon, urcrnrlon=maxlon, resolution='c', ax=self.axes)

#        basemap = Basemap(projection='tmerc', llcrnrlat=minlat, urcrnrlat=maxlat,\
#                    llcrnrlon=minlon, urcrnrlon=maxlon, resolution='c', lon_0=-4.36,
#                     lat_0=54.7,ax=self.axes)
#        basemap = Basemap(projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,\
#                    llcrnrlon=minlon, urcrnrlon=maxlon, resolution='c', ax=self.axes)
#        basemap.drawcoastlines()
        # draw parallels and meridians.
#        basemap.drawparallels(np.arange(-90., 91., 30.))
#        basemap.drawmeridians(np.arange(-180., 181., 60.))
#        basemap.drawmapboundary(fill_color='aqua')
        clevs = [-1, 0, 1]
        bathy_clevs = np.arange(-1, 5300, 100)
        x = np.arange(0, self.mask.lon.shape[0])
        y = np.arange(0, self.mask.lon.shape[1])
        x_vals, y_vals = np.meshgrid(y, x)
#        x_vals = self.mask.lon
#        y_vals = self.mask.lat
        Z = self.mask.bathy_data[...].astype(np.float64)
        self.axes.contourf(x_vals, y_vals, Z, 10, cmap=plt.get_cmap('GnBu'))#cmap=cm.s3pcpn)
        self.axes.contourf(x_vals, y_vals, self.mask.data, 5, cmap=plt.get_cmap('autumn'),\
                           alpha=mask_alpha)
#        basemap.contourf(x_vals, y_vals, )
#        plt.title("Equidistant Cylindrical Projection")
        self.canvas.draw()

    def add_mask(self):
        """ adds the selected region in the drawing tool to the mask """
        if self._drawing_tool_name != "" and self.mask != None:
            if self._drawing_tool.polygon != None:
#                grid = zip(self.mask.lon.ravel(), self.mask.lat.ravel())
                x = np.arange(0, self.mask.lon.shape[0])
                y = np.arange(0, self.mask.lon.shape[1])
                x_vals, y_vals = np.meshgrid(y, x)
                grid = zip(x_vals.ravel(), y_vals.ravel())

                self._drawing_tool.polygon.set_linewidth(1.0)
                p_path = Path(self._drawing_tool.polygon.xy)
                index = p_path.contains_points(grid)
                index = index.reshape(self.mask.lon.shape)
                tmp = self.mask.data[index]
                tmp[tmp == 1] = -1
                self.mask.data[index] = tmp
                self._drawing_tool.reset()
                self.axes.clear()
                self.create_basemap()

    def remove_mask(self):
        """ removes the selected region in the drawing tool from the mask """
        if self._drawing_tool_name != "" and self.mask != None:
            if self._drawing_tool.polygon != None:
#                grid = zip(self.mask.lon.ravel(), self.mask.lat.ravel())
                x = np.arange(0, self.mask.lon.shape[0])
                y = np.arange(0, self.mask.lon.shape[1])
                x_vals, y_vals = np.meshgrid(y, x)
                grid = zip(x_vals.ravel(), y_vals.ravel()) #check for the index

                self._drawing_tool.polygon.set_linewidth(1.0)
                p_path = Path(self._drawing_tool.polygon.xy)
                index = p_path.contains_points(grid)
                index = index.reshape(self.mask.lon.shape)
                tmp = self.mask.data[index]
                tmp[tmp == -1] = 1
                self.mask.data[index] = tmp
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

    @pyqtSlot(str,str)
    def set_bathymetry_file(self, bathymetry_filename, mask_file):
        """ Set the bathymetry file """
        self.mask = Mask(bathymetry_filename, mask_file)
        self.create_basemap()
        
    @pyqtSlot(str)
    def save_mask_file(self, mask_file):
        """ Save the mask data to mask_file """
        if self.mask is not None:
            self.mask.save_mask(mask_file)

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
                          (None, None, None, None),\
                          ('Freehand', 'Freehand drawing', 'freehand', 'freehand'),\
                          ('Rectangle', 'Rectangle drawing', 'rectangle', 'rectangle'),\
                          ('Border', 'Border selection', 'border', 'border'),\
                          ('plus', 'Add mask', 'add_mask', 'add_mask'),\
                          ('minus', 'Remove mask', 'remove_mask', 'remove_mask'),\
                          (None, None, None, None)\
                          )
        NavigationToolbar.__init__(self, canvas, parent)
        self._actions['freehand'].setCheckable(True)
        self._actions['freehand'].setIcon(set_icon('freehand.png'))
        self._actions['rectangle'].setCheckable(True)
        self._actions['rectangle'].setIcon(set_icon('rectangle.png'))
        self._actions['border'].setIcon(set_icon('border.png'))
        self._actions['add_mask'].setIcon(set_icon('plus.png'))
        self._actions['remove_mask'].setIcon(set_icon('minus.png'))


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

def set_icon(name):
    """ Creates an icon based on the file found in the module directory with input name"""
    return QtGui.QIcon(os.path.join(os.path.dirname(__file__), name))

#if __name__ == '__main__':
#    mask = Mask('..\..\data\grid_C\NNA_R12_bathy_meter_bench.nc')
#    app = QtGui.QApplication(sys.argv)
#    widget = MatplotlibWidget()
#    widget.set_bathymetry_file('..\..\data\grid_C\NNA_R12_bathy_meter_bench.nc')
#    #mask=mask)
#    fig = widget.figure
#    ax = widget.axes#fig.add_subplot(111)
#     map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
#             llcrnrlon=-180,urcrnrlon=180,resolution='c', ax=ax)
#     map.drawcoastlines()
#     map.fillcontinents(color='coral',lake_color='aqua')
#         # draw parallels and meridians.
#     map.drawparallels(np.arange(-90.,91.,30.))
#     map.drawmeridians(np.arange(-180.,181.,60.))
#     map.drawmapboundary(fill_color='aqua')
#    fig.canvas.draw()
#    widget.show()
#    app.exec_()
