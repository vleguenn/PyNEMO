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
   sn_src_hgr = 'F:/NEMO_bdy_tools/bdy_matlab/grid_low_res_C/mesh_hgr.nc'   !  /grid/
   sn_src_zgr = 'F:/NEMO_bdy_tools/bdy_matlab/grid_low_res_C/mesh_zgr.nc'
   sn_dst_hgr = 'F:/NEMO_bdy_tools/bdy_matlab/grid_C/mesh_hgr_zps.nc'
   sn_dst_zgr = 'F:/NEMO_bdy_tools/bdy_matlab/grid_C/mesh_zgr_zps.nc'
   sn_src_msk = 'F:/NEMO_bdy_tools/bdy_matlab/grid_low_res_C/mask.nc'

!-----------------------------------------------------------------------
!  I/O 
!-----------------------------------------------------------------------
   sn_src_dir = 'F:/NEMO_bdy_tools/bdy_matlab/srcdata_low_res_C/'       ! src_files/'
   sn_dst_dir = 'F:/NEMO_bdy_tools/scratch/jofa/'
   sn_fn      = 'NNA_R12'                 ! prefix for output files
   nn_fv      = -1e20                     !  set fill value for output files
   nn_src_time_adj = 0					  ! src time adjustment

!-----------------------------------------------------------------------
!  unstructured open boundaries                         
!-----------------------------------------------------------------------
    ln_coords_file = .true.               !  =T : produce bdy coordinates files
    cn_coords_file = 'coordinates.bdy.nc' !  name of bdy coordinates files (if ln_coords_file=.TRUE.)
    ln_mask_file   = .false.              !  =T : read mask from file
    cn_mask_file   = ''                   !  name of mask file (if ln_mask_file=.TRUE.)
    ln_dyn2d       = .true.               !  boundary conditions for barotropic fields
    ln_dyn3d       = .true.               !  boundary conditions for baroclinic velocities
    ln_tra         = .true.               !  boundary conditions for T and S
    ln_ice         = .false.               !  ice boundary condition   
    nn_rimwidth    = 9                    !  width of the relaxation zone

!-----------------------------------------------------------------------
!  Time information
!-----------------------------------------------------------------------
    nn_year_000     = 1982        !  year start
    nn_year_end     = 1982        !  year end
    nn_month_000    = 1           !  month start (default = 1 is years>1)
    nn_month_end    = 12          !  month end (default = 12 is years>1)
    sn_dst_calendar = 'gregorian' !  output calendar format
    nn_base_year    = 1960        !  base year for time counter

!-----------------------------------------------------------------------
!  Additional parameters
!-----------------------------------------------------------------------
    nn_wei  = 1                   !  smoothing filter weights 
    rn_r0   = 0.0417              !  decorrelation distance use in gauss
                                  !  smoothing onto dst points. Need to 
                                  !  make this a funct. of dlon
    sn_history  = 'EB bdy files produced by jofa from ORCA0083-N001 for testing'
                                  !  history for netcdf file
    ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
    nn_alpha    = 0               !  Euler rotation angle
    nn_beta     = 0               !  Euler rotation angle
    nn_gamma    = 0               !  Euler rotation angle
