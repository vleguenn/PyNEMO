'''
funtion e3_to_depth
Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp.
Action  :   pe3t, pe3w : scale factor of t- and w-point (m)
Useage: [gdept, gdepw] = e3_to_depth(e3t, e3w, nz)
'''

import numpy as np

def e3_to_depth(pe3t, pe3w, jpk):
  pdepw      = np.zeros_like(pe3w)
  pdepw[0,:] = 0.
  pdept      = np.zeros_like(pe3t)
  pdept[0,:] = 0.5 * pe3w[0,:]

  for jk in np.arange(1,jpk,1):
    pdepw[jk,:] = pdepw[jk-1,:] + pe3t[jk-1,:]
    pdept[jk,:] = pdept[jk-1,:] + pe3w[jk  ,:]

  return pdept, pdepw
