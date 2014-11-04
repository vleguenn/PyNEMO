############################################################
# Nemo Boundary Generation                                 #
# Creates indices for t, u and v points, plus rim gradient #
# # # #                                              # # # #
# Written by John Kazimierz Farey in August 2012           #
# Port of Matlab code by James Harle                       #
############################################################


"""
I have renamed some equivalent variables from matlab to try to follow 
the python style guide- Also bdy_i = bdy_t: I thought it added clarity
since this class will generates the u and v points as well as t. 

Initial arguments:
boundary_mask: array
settings: dictionary
grid: string 't', 'u', 'v' or 'f' specifies grid type

Attributes:
Initially: bdy_i, bdy_r
Additional assigned externally

"""

import numpy as np
import logging

from nemo_bdy_lib import sub2ind

class Boundary:

    # Bearings for overlays
    _NORTH = [1,-1,1,-1,2,None,1,-1]
    _SOUTH = [1,-1,1,-1,None,-2,1,-1]
    _EAST = [1,-1,1,-1,1,-1,2,None]
    _WEST = [1,-1,1,-1,1,-1,None,-2]
    
    def __init__(self, boundary_mask, settings, grid):
        self.logger = logging.getLogger(__name__)
        bdy_msk = boundary_mask
        self.settings = settings
        rw = self.settings['rimwidth']
        rm = rw - 1
        self.grid_type = grid.lower()
        # Throw an error for grid type confusion
        if grid not in ['t', 'u', 'v', 'f']:
            self.logger.error('Grid Type not correctly specified:'+grid)
            raise ValueError("""%s is invalid grid type;  
                                must be 't', 'u', 'v' or 'f'""" %grid)

        # Configure grid type
        if grid is not 't':
            # We need to copy this before changing, because the original will be 
            # needed in calculating later grid boundary types
            bdy_msk = boundary_mask.copy()
            grid_ind = np.zeros(bdy_msk.shape, dtype=np.bool, order='C')
            for fval in [-1, 0]:
                if grid is 'u':
                    grid_ind[:, :-1] = np.logical_and(bdy_msk[:, :-1] == 1, 
                                                       bdy_msk[:, 1:] == fval)
                    bdy_msk[grid_ind] = fval
                elif grid is 'v':
                    grid_ind[:-1, :] = np.logical_and(bdy_msk[:-1, :] == 1, 
                                                       bdy_msk[1:, :] == fval)
                    bdy_msk[grid_ind] = fval
                elif grid is 'f':
                    grid_ind[:-1, :-1] = np.logical_and(bdy_msk[:-1,:-1] == 1,
                                                       bdy_msk[1:, 1:] == fval)
                    grid_ind[:-1, :] = np.logical_or(
                                        np.logical_and(bdy_msk[:-1, :] == 1,
                                                       bdy_msk[1:, :] == fval), 
                                                     grid_ind[:-1, :] == 1)

                    grid_ind[:, :-1] = np.logical_or(
                                        np.logical_and(bdy_msk[:, :-1] == 1,
                                                       bdy_msk[:, 1:] == fval), 
                                                     grid_ind[:, :-1] == 1)
                    bdy_msk[grid_ind] = fval

        # Create padded array for overlays
        msk = self._pad_2d(bdy_msk, -1)
        
        # create index arrays of I and J coords
        igrid, jgrid = np.meshgrid((np.arange(len(bdy_msk[0,:]))), 
                                   (np.arange(len(bdy_msk[:,0]))))

        SBi, SBj = self._find_bdy(igrid, jgrid, msk, self._SOUTH)
        NBi, NBj = self._find_bdy(igrid, jgrid, msk, self._NORTH)
        EBi, EBj = self._find_bdy(igrid, jgrid, msk, self._EAST)
        WBi, WBj = self._find_bdy(igrid, jgrid, msk, self._WEST)

        ti = np.r_[SBi,NBi,WBi,EBi]
        tj = np.r_[SBj,NBj,WBj,EBj]
        tij = np.c_[ti,tj]
        bdy_i = np.tile(tij, (rw,1,1))
        
        bdy_i = np.transpose(bdy_i, (1,2,0))
        bdy_r = np.squeeze(bdy_i[:,0,:]) * 0

        # Add points for relaxation zone over rim width
        temp = np.array([np.r_[SBi*0, NBi*0, WBi*0+1, EBi*0-1], 
                        np.r_[SBj*0+1, NBj*0-1, WBj*0, EBj*0]])
        for i in range(rm):
            bdy_i[:,:,i+1] = bdy_i[:,:,i] + np.transpose(temp, (1,0))
            bdy_r[:,i+1] = i+1
                   
        bdy_i = np.transpose(bdy_i, (1,2,0))
        bdy_i = np.reshape(bdy_i, 
                 (bdy_i.shape[0],bdy_i.shape[1]*bdy_i.shape[2]))
        bdy_r = bdy_r.flatten(1)

        ##   Remove duplicate and open sea points  ##

        # Note: Needs to use a stable sort method to preserve order
        # Get index of ascending order. 
        s_ind = np.argsort(bdy_r,kind='mergesort')

        # Rearrange according to the order
        bdy_r = bdy_r[s_ind]
        bdy_i = bdy_i[:,s_ind]
        
        # copy might not be ideal
        bdy_i2 = np.transpose(bdy_i,(1,0)).copy()
        uniqind = self._unique_rows(bdy_i2)
        bdy_i = bdy_i2[uniqind]
        bdy_r = bdy_r[uniqind]
        
        ###   Fill in any gradients between relaxation zone and internal domain
        ###   bdy_msk matches matlabs incarnation, r_msk is pythonic 
        r_msk = bdy_msk.copy()
        r_msk[r_msk == 1] = rw
        r_msk = np.float16(r_msk)
        r_msk[r_msk < 1] = np.NaN
        
        grad_ind = sub2ind((bdy_msk.shape[0], bdy_msk.shape[1]), 
                    bdy_i[:,0], bdy_i[:,1])

        lm_ind = bdy_msk.flatten(1)[grad_ind] == 0
        bdy_i = bdy_i[(np.invert(lm_ind)),:]
        bdy_r = bdy_r[np.invert(lm_ind)]
        grad_ind = grad_ind[np.invert(lm_ind)]
        
        rm_shape = r_msk.shape
        r_msk = r_msk.flatten(1)
        r_msk[grad_ind] = np.float16(bdy_r)
        r_msk = r_msk.reshape(rm_shape, order='F')

        r_msk_orig = r_msk.copy()
        r_msk_ref = r_msk[1:-1,1:-1]

        self.logger.debug('Start r_msk bearings loop')
        # # # This is by far the slowest part # # #
        for i in range(rw-1):
            # Check each bearing
            for b in [self._SOUTH, self._NORTH, self._WEST, self._EAST]:
                r_msk,r_msk_ref = self._fill(r_msk,r_msk_ref, b)
        self.logger.debug('done loop')
    
        # update bdy_i and bdy_r
        new_ind = np.abs(r_msk - r_msk_orig) >  0
        # The transposing gets around the Fortran v C ordering thing.
        bdy_i_tmp = np.array([igrid.T[new_ind.T], jgrid.T[new_ind.T]])
        bdy_r_tmp = r_msk.T[new_ind.T]
        bdy_i = np.vstack((bdy_i_tmp.T, bdy_i))
        
        uniqind = self._unique_rows(bdy_i)
        bdy_i = bdy_i[uniqind, :]
        bdy_r = np.hstack((bdy_r_tmp, bdy_r))
        bdy_r = bdy_r[uniqind]
        
        # sort by rimwidth
        igrid = np.argsort(bdy_r,kind='mergesort')
        bdy_r = bdy_r[igrid]
        bdy_i = bdy_i[igrid,:]

        self.bdy_i = bdy_i
        self.bdy_r = bdy_r

        self.logger.debug( 'Final bdy_i: %s', self.bdy_i.shape)


# # # # # # # # # # # # #
# # # # Functions # # # # 
# # # # # # # # # # # # # 

    # pads each edge of 2d array with val
    def _pad_2d(self, arr, val):
        arr = np.insert(arr,0,val,0)
        arr = np.insert(arr,0,val,1)
        arr = np.insert(arr,arr.shape[0],val,0)
        arr = np.insert(arr,arr.shape[1],val,1)

        return arr

    def _find_bdy(self, I,J, mask, brg):
        # subtract matrices to find boundaries, set to True
        m1 = mask[brg[0]:brg[1], brg[2]:brg[3]]
        m2 = mask[brg[4]:brg[5], brg[6]:brg[7]]
        overlay = np.subtract(m1,m2)
        # Create boolean array of bdy points in overlay
        bool_arr = overlay==2
        # index I or J to find bdies
        bdy_I = I[bool_arr]
        bdy_J = J[bool_arr]
        
        return bdy_I, bdy_J

    # Find unique rows in an array
    def _unique_rows(self, arr):
        a_str = np.empty([len(arr[:,0]), 1], dtype='|S8')
        for ii in range(len(arr[:,0])):
            a_str[ii] = '%s' %(''.join('%4.4i' % i for i in arr[ii,:]))
        uniq, ind =  np.unique(a_str, return_index=True)

        return  ind

    def _fill(self, mask, ref, brg):
        tmp = mask[brg[4]:brg[5],brg[6]:brg[7]]
        ind = (ref - tmp) > 1
        ref[ind] = tmp[ind] + 1
        mask[brg[0]:brg[1],brg[2]:brg[3]] = ref

        return mask, ref


    """
    # Need Numpy 1.7 mergesort support to return index on unique func- 
    # v1.6x does not support all dtypes. If we could use this it might
    # be much faster
    def alt_unique_rows(t):
        
        tlist = t.tolist()
        sortt = []
        indx = zip(*sorted([(val, i) for i,val in enumerate(tlist)]))[1]
        indx = np.array(indx)
        for i in indx:
            sortt.append(tlist[i])
        del tlist
        for i,x in enumerate(sortt):
            if x == sortt[i-1]:
                indx[i] = -1
        
        return indx[indx != -1]

    def sub2ind(self, shap, subx, suby):
        ind = (subx * shap[0]) + suby
        #ind = (suby % shap[0]) * shap[1] + (subx % shap[1])
        return ind

    """


