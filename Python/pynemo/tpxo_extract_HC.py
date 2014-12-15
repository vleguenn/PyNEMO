'''
This is to extract the tidal harmonic constants out of a tidal model
for a given locations
[amp,Gph] = tpxo_extract_HC(Model,lat,lon,type,Cid)

@author: Mr. Srikanth Nagella
'''

from netCDF4 import Dataset
from scipy import interpolate
import numpy as np
from astropy.units import amp
from lib2to3.pgen2.token import AMPER
 
class TPXO_Extract:
    def __init__(self,lat,lon,type):
        #read the grid file
        self.grid = Dataset('../data/tide/grid_tpxo7.2.nc')
        #read the h file
        self.h = Dataset('../data/tide/h_tpxo7.2.nc')
        #read the u file
        self.u = Dataset('../data/tide/u_tpxo7.2.nc')

        H = self.grid.variables['hz']
        mz = self.grid.variables['mz']
        x  = self.h.variables['lon_z'][:,0]
        y  = self.h.variables['lat_z'][0,:]
        dx = x[1] - x[0] 
        dy = y[1] - y[0]
        km = 0 # added to maintain the reference to matlab tmd code    
        # Pull out the constituents that are avaibable
        self.cons = []
        for ncon in range(self.h.variables['con'].shape[0]):
            self.cons.append(self.h.variables['con'][ncon,:].tostring().strip())
        print self.cons

        glob = 0
        if x[-1]-x[0] == 360-dx:
            glob = 1
        
        if glob == 1:
            x = np.concatenate(([x[0]-dx,],x,[x[-1]+dx,]))
            H = np.concatenate(([H[-1,:],],H,[H[0,:],]),axis=0)
            mz = np.concatenate(([mz[-1,:],],mz,[mz[0,:],]),axis=0)
    
        #adjust lon convention
        xmin = np.min(lon)
        xmax = np.max(lon)
        if km ==0 :
            if xmin < x[0]:
                lon[lon<0] = lon[lon<0] + 360
            if xmin > x[-1]:
                lon[lon>180] = lon[lon>180]-360

        #H[H==0] = np.NaN
#        f=interpolate.RectBivariateSpline(x,y,H,kx=1,ky=1)
#        D = np.zeros(lon.size) 
#        for idx in range(lon.size):
#            D[idx] = f(lon[idx],lat[idx])
#        print D[369:371]

#        H2 = np.ravel(H)
#        H2[H2==0] = np.NaN
#        points= np.concatenate((np.ravel(self.h.variables['lon_z']),np.ravel(self.h.variables['lat_z'])))
#        points= np.reshape(points,(points.shape[0]/2,2),order='F')
#        print points.shape
#        print np.ravel(H).shape
#        D = interpolate.griddata(points,H2,(lon,lat))
#        print D
#        print D.shape

        H[H==0] = np.NaN
        lonlat = np.concatenate((lon,lat))
        lonlat = np.reshape(lonlat,(lon.size,2),order='F')

        D = interpolate.interpn((x,y), H, lonlat)
     
#        f=interpolate.RectBivariateSpline(x,y,mz,kx=1,ky=1)
#        mz1 = np.zeros(lon.size)
#        for idx in range(lon.size):
#            mz1[idx] = f(lon[idx],lat[idx])
        mz1 = interpolate.interpn((x,y),mz,lonlat)        
        
        i1 = np.where((np.isnan(D)) & (mz1 > 0))
        if i1[0].size!=0:
            D[i1] = self.BLinterp(x,y,H, lon[i1], lat[i1])
        
        if type == 'z' or type== 't':
            self.amp, self.Gph = self.interpolate_constituents(self.h, 'hRe', 'hIm', 'lon_z','lat_z', lon, lat,maskname='mz')
        elif type == 'u':
            self.amp, self.Gph = self.interpolate_constituents(self.u, 'URe', 'UIm', 'lon_u', 'lat_u', lon, lat,D,maskname='mu')                              
        elif type == 'v':
            self.amp, self.Gph = self.interpolate_constituents(self.u, 'VRe', 'VIm', 'lon_v', 'lat_v', lon, lat,D,maskname='mv')
        else:
            print 'Unknown type'
            return    
        print self.amp
        print self.Gph
    
    def interpolate_constituents(self, nc, realName, imgName, lonName, latName, lon, lat, D=None, maskname=None):
        amp = np.zeros((nc.variables['con'].shape[0], lon.shape[0],))
        Gph = np.zeros((nc.variables['con'].shape[0], lon.shape[0],))            
        z = np.array(np.ravel(nc.variables[realName]), dtype=complex)
        z.imag = np.array(np.ravel(nc.variables[imgName]))
        z = z.reshape(nc.variables[realName].shape)
        
        print z.shape
        #z[z==0] = np.NaN
        
        #Lat Lon values
        x = nc.variables[lonName][:,0]
        y = nc.variables[latName][0,:]
        dx = x[1] - x[0] 
        dy = y[1] - y[0]
        glob = 0
        if x[-1]-x[0] == 360-dx:
            glob = 1
        
        if glob == 1:
            x = np.concatenate(([x[0]-dx,],x,[x[-1]+dx,]))
    
        #adjust lon convention
        xmin = np.min(lon)
        xmax = np.max(lon)
        
        if xmin < x[0]:
            lon[lon<0] = lon[lon<0] + 360
        if xmin > x[-1]:
            lon[lon>180] = lon[lon>180]-360

        lonlat = np.concatenate((lon,lat))
        lonlat = np.reshape(lonlat,(lon.size,2),order='F')
        
                       
        mask = self.grid.variables[maskname]
        mask = np.concatenate(([mask[-1,:],],mask,[mask[0,:],]),axis=0)
        #interpolate the mask values
        maskedpoints = interpolate.interpn((x,y),mask,lonlat)

        z1 = np.zeros((z.shape[0], lon.shape[0], 2,))
        for ic in range(z.shape[0]):
            #interpolate real values
            z1[ic,:,0] = self.interpolate_data(x, y, z[ic,:,:].real, maskedpoints, lonlat)
            #interpolate imag values
            z1[ic,:,1] = self.interpolate_data(x, y, z[ic,:,:].imag, maskedpoints, lonlat)
            
            #for velocity values
            if D is not None:
                z1[ic,:,0] = z1[ic,:,0]/D*100
                z1[ic,:,1] = z1[ic,:,1]/D*100

            zcomplex = np.array(z1[ic, :, 0], dtype=complex)
            zcomplex.imag = z1[ic, :, 1]

            amp[ic, :] = np.absolute(zcomplex)
            Gph[ic, :] = np.arctan2(-1*zcomplex.imag, zcomplex.real)
        Gph=Gph*180.0/np.pi
        Gph[Gph<0] = Gph[Gph<0]+360.0
        return amp, Gph
            
    def interpolate_data(self,x,y,H,m,lonlat):
        z1 = np.zeros((lonlat.shape[0],))
        H[H==0]=np.NaN
        H = np.concatenate(([H[-1,:],],H,[H[0,:],]),axis=0)
        z1[:] = interpolate.interpn((x,y), H, lonlat)            
        i1 = np.where((np.isnan(z1)) & (m > 0))
        if i1[0].size!=0:
            z1[i1] = self.BLinterp(x,y,H, np.ravel(lonlat[i1,0]), np.ravel(lonlat[i1,1]))
        return z1
                    
    def BLinterp(self,x,y,h,xt,yt):
        glob = 0
        dx=x[1]-x[0]
        if x[-1] - x[1]==360-dx:
            glob = 1
        inan = np.where(np.isnan(h))
        h[inan] = 0
        mz = np.zeros(h.shape)
        mz[h!=0]=1
        n = x.size
        m = y.size
        if n!=h.shape[0] or m!=h.shape[1]:
            print 'Check Dimensions'
            return np.NaN
        if glob == 1:
            h0=h
            mz0=mz
            x=np.concatenate(([x[0]-2*dx,x[0]-dx,],x,[x[-1]+dx,x[-1]+2*dx]))
            h=np.concatenate((h[-2,:],h[-1,:],h,h[0,:],h[1,:]),axis=0)
            mz=np.concatenate((mz[-2,:],mz[-1,:],mz,mz[0,:],mz[1,:]),axis=0)
        xti = xt
        
        ik = np.where((xti<x[0]) & (xti>x[-1]))
        if x[-1] > 180:
            xti[ik] = xti[ik]+360
        if x[-1] < 0:
            xti[ik] = xti[ik]-360
        xti[xti>360] = xti[xti>360]-360
        xti[xti<-180] = xti[xti<-180]+360        
        
        q=1/(4+2*np.sqrt(2))
        q1 = q/np.sqrt(2)
        h1 = q1*h[0:-2,0:-2]+q*h[0:-2,1:-1]+q1*h[0:-2,2:]+q1*h[2:,0:-2]+q*h[2:,1:-1]+q1*h[2:,2:]+q*h[1:-1,0:-2]+q*h[1:-1,2:]
        mz1 = q1*mz[0:-2,0:-2]+q*mz[0:-2,1:-1]+q1*mz[0:-2,2:]+q1*mz[2:,0:-2]+q*mz[2:,1:-1]+q1*mz[2:,2:]+q*mz[1:-1,0:-2]+q*mz[1:-1,2:]
        mz1[mz1==0]=1
        h2 = h.copy()
        h2[1:-1,1:-1] = np.divide(h1,mz1)
        ik = np.where(mz==1)
        h2[ik] = h[ik]
        h2[h2==0] = np.NaN
        lonlat = np.concatenate((xti,yt))
        lonlat = np.reshape(lonlat,(xti.size,2),order='F')        
        hi = interpolate.interpn((x,y),h2,lonlat)
        
        return hi
          
#lat=[42.8920,42.9549,43.0178]
#lon=[339.4313,339.4324,339.4335]
#lat_u=[42.8916,42.9545,43.0174]
#lon_u=[339.4735,339.4746,339.4757]
#lat = np.array(lat_u)
#lon = np.array(lon_u)
#x = TPXO_Extract(lat,lon,'u')
