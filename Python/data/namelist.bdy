!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OPA  : namelist for BDY generation tool
!!            
!!             User inputs for generating open boundary conditions
!!             employed by the BDY module in NEMO. Boundary data
!!             can be set up for v3.2 NEMO and above.
!!            
!!             More info here.....
!!            
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!-----------------------------------------------------------------------
!   vertical coordinate
!-----------------------------------------------------------------------
   ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)  
   ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
   ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
   rn_hmin     =   -10     !  min depth of the ocean (>0) or 
                           !  min number of ocean level (<0)

!-----------------------------------------------------------------------
!   s-coordinate or hybrid z-s-coordinate
!-----------------------------------------------------------------------
   rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
   rn_sbot_max = 7000.     !  maximum depth of s-bottom surface 
                           !  (= ocean depth) (>0) (m)
   ln_s_sigma  = .false.   !  hybrid s-sigma coordinates
   rn_hc       =  150.0    !  critical depth with s-sigma

!-----------------------------------------------------------------------
!  grid information 
!-----------------------------------------------------------------------
   sn_src_hgr = 'http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/grid_low_res_C/mesh_hgr.nc'   !  /grid/
   sn_src_zgr = 'http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/grid_low_res_C/mesh_zgr.nc'
   sn_dst_hgr = 'http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/grid_C/mesh_hgr_zps.nc'
   sn_dst_zgr = 'http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/grid_C/mesh_zgr_zps.nc'
   sn_src_msk = 'http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/grid_low_res_C/mask.nc'
   sn_bathy   = 'http://esurgeod.noc.soton.ac.uk:8080/thredds/dodsC/PyNEMO/grid_C/NNA_R12_bathy_meter_bench.nc'

!-----------------------------------------------------------------------
!  I/O 
!-----------------------------------------------------------------------
!   sn_src_dir = 'F:/NEMO_bdy_tools/bdy_matlab/srcdata_low_res_C/'       ! src_files/'
!   sn_src_dir = 'http://esurgeod.noc.soton.ac.uk:8080/thredds/catalog/PyNEMO/data/catalog.xml'
   sn_src_dir = '/Users/srikanthnagella/Development/Python/pynemo/pynemo/trunk/Python/data/test2.ncml'       ! src_files/'
   sn_dst_dir = 'scratch/jofa/'
   sn_fn      = 'NNA_R12'                 ! prefix for output files
   nn_fv      = -1e20                     !  set fill value for output files
   nn_src_time_adj = 0					  ! src time adjustment
   sn_dst_metainfo = 'EB bdy files produced by jdha from ORCA0083-N001 gloabl run provided by acc'

!-----------------------------------------------------------------------
!  unstructured open boundaries                         
!-----------------------------------------------------------------------
    ln_coords_file = .true.               !  =T : produce bdy coordinates files
    cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
    ln_mask_file   = .false.              !  =T : read mask from file
    cn_mask_file   = 'F:/NEMO_bdy_tools/scratch/testmask.nc'                   !  name of mask file (if ln_mask_file=.TRUE.)
    ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
    ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
    ln_tra         = .true.               !  boundary conditions for T and S
    ln_ice         = .true.               !  ice boundary condition   
    nn_rimwidth    = 9                    !  width of the relaxation zone

!-----------------------------------------------------------------------
!  unstructured open boundaries tidal parameters                        
!-----------------------------------------------------------------------
    ln_tide        = .false.               !  =T : produce bdy tidal conditions
    clname(1)      = 'M2'                 ! constituent name
    clname(2)      = 'S2'         
    clname(3)      = 'K2'        
    ln_trans       = .true.

!-----------------------------------------------------------------------
!  Time information
!-----------------------------------------------------------------------
    nn_year_000     = 1979        !  year start
    nn_year_end     = 1979        !  year end
    nn_month_000    = 11           !  month start (default = 1 is years>1)
    nn_month_end    = 11          !  month end (default = 12 is years>1)
    sn_dst_calendar = 'gregorian' !  output calendar format
    nn_base_year    = 1960        !  base year for time counter
	sn_tide_grid   = 'F:/NEMO_bdy_tools/bdy_matlab/bdy_matlab/tmd/DATA/grid_tpxo7.2.nc'
	sn_tide_h	   = 'F:/NEMO_bdy_tools/bdy_matlab/bdy_matlab/tmd/DATA/h_tpxo7.2.nc'
	sn_tide_u	   = 'F:/NEMO_bdy_tools/bdy_matlab/bdy_matlab/tmd/DATA/u_tpxo7.2.nc'
	
!-----------------------------------------------------------------------
!  Additional parameters
!-----------------------------------------------------------------------
    nn_wei  = 1                   !  smoothing filter weights 
    rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                  !  smoothing onto dst points. Need to 
                                  !  make this a funct. of dlon
    sn_history  = 'EB bdy files produced by jofa from ORCA0083-N001 for testing'
                                  !  history for netcdf file
    ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
    nn_alpha    = 0               !  Euler rotation angle
    nn_beta     = 0               !  Euler rotation angle
    nn_gamma    = 0               !  Euler rotation angle
	rn_mask_max_depth = 100.0	  !  Maximum depth to be ignored for the mask
	rn_mask_shelfbreak_dist = 20000.0    !  Distance from the shelf break
