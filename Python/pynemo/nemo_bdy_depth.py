#
#
#
#
#
#


class Depth:

    def __init__(self):
        def gen_depth(self), bdy_ind_i, bdy_ind_j, hmin): #, sco, hc):
        '''
        #dst_zgr = '/work/jofa/grid/mesh_zgr_zps.nc'
        bdy_ind_i = self.bdy_i[:,0]
        bdy_ind_j = self.bdy_i[:,1]
        nbdy = len(bdy_ind_i)
        '''


        sco = self.settings['sco'] #False
        hc = self.settings['hc']#150

        nc = Dataset(self.settings['dst_zgr'], 'r')
        mbathy = nc.variables['mbathy'][:,:,:]
        # numpy can't NaN ints. question this
        mbathy = np.float16(mbathy)
        mbathy[mbathy is 0] = np.NaN 
        nz = len(nc.variables['nav_lev'][:])
        

        '''
        # Set up arrays
        self.z_points = np.zeros((nz, nbdy))
        self.z_w_points = np.zeros((nz, nbdy))
        '''     
        # Check inputs
        # FIX ME? Errors for wrong obj arg len. probably better to work around
        if sco:
            # hc = ... FIX ME??
            # Depth of water column at t-point
            hbatt = nc['hbatt'][:,:,:]
            # Replace land with NaN   
            hbatt[mbathy is 0] = np.NaN
        '''
        # find bdy indices from subscripts
        ind = self.sub2ind(mbathy.shape, bdy_ind_j, bdy_ind_i)
        ind2 = 0
        if self.grid_type is 'u':
            ind2 = self.sub2ind(mbathy.shape, bdy_ind_j, bdy_ind_i + 1)
        elif self.grid_type is 'v':
            ind2 = self.sub2ind(mbathy.shape, bdy_ind_j + 1, bdy_ind_i)
        '''
        for k in range(nz):
            if sco:
                # sigma coeffs at t-point (1->0 indexed)
                gsigt = nc['gsigt'][0,k,:,:]
                # sigma coeffs at w-point
                gsigw = nc['gsigw'][0,k,:,:]

                # check size of gsigt SKIPPED

                wrk1 = (hbatt - hc) * gsigt[:,:] + (hc * (k - 0.5) / (nz - 1))
                wrk2 = (hbatt - hc) * gsigw[:,:] + (hc * (k - 0.5) / (nz - 1))
            else:
                wrk1 = nc['gdept'][0,k,:,:]
                wrk2 = nc['gdepw'][0,k,:,:]

            # Replace deep levels that are not used with NaN
            wrk2[mbathy + 1 < k] = np.NaN
`           wrk1[mbathy < k] = np.NaN


            # Set u and v grid point depths
            if ind2:
                self.z_w_points[k,:] = 0.5 * (wrk2[ind] + wrk2[ind2]) ; 
                self.z_points[k,:]  = 0.5 * (wrk1[ind] + wrk1[ind2]) ; 
            
            else:           
                self.z_w_points[k,:] = wrk2[ind] ;
                self.z_points[k,:]  = wrk1[ind] ;        
                                

        nc.close()



