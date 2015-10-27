Examples
========
Here we provide two worked examples using pyNEMO. The first is a setup of the Northwest European Shelf using both local
and remote dataset. The second is an end-to-end setup of a small regional model in the tropics.

Example 1: Northwest European Shelf
-----------------------------------

Example 2: Lighthouse Reef
--------------------------

This example has been tested on the ARCHER HPC facillity.
First, create a working directory into which the NEMO 
source code can be checked out. Create an inputs directory
to unpack the forcing tar ball.

::

   cd $WDIR
   mkdir INPUTS
   cd INPUTS
   wget INPUTS.tar.gz
   tar xvfz INPUTS.tar.gz
   rm INPUTS.tar.gz
   cd ../
   svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/2014/dev_r4621_NOC4_BDY_VERT_INTERP


Next we setup our experiment directory (we are only using 
OPA_SRC) and drop an updated dtatsd.F90 into MY_SRC to 
allow the vertical interpolation of initial conditions on 
to the new verictal coordinates. We also apply a patch to 
fldread.F90 for missing variables.

::

   export CDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/CONFIG
   export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS
   cd $CDIR/../NEMO/OPA_SRC/SBC
   patch -b < $WDIR/fldread.patch
   cd $CDIR
   ./makenemo -n LH_REEF -m XC_ARCHER_INTEL
   cp $WDIR/INPUTS/dtatsd.F90 LH_REEF/MY_SRC/ 

To generate bathymetry, initial conditions and grid information
we first need to compile some of the NEMO TOOLS (after a small
bugfix - and to allow direct passing of arguments). For some 
reason GRIDGEN doesn't like INTEL:

::

   cd $WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS/WEIGHTS/src
   patch -b < $WDIR/INPUTS/scripinterp_mod.patch
   patch -b < $WDIR/INPUTS/scripinterp.patch
   patch -b < $WDIR/INPUTS/scrip.patch
   patch -b < $WDIR/INPUTS/scripshape.patch
   patch -b < $WDIR/INPUTS/scripgrid.patch
   cd ../../
   cp $WDIR/INPUTS/arch-* ../ARCH
   ./maketools -n WEIGHTS -m XC_ARCHER_INTEL
   ./maketools -n REBUILD_NEMO -m XC_ARCHER_INTEL
   module unload cray-netcdf-hdf5parallel cray-hdf5-parallel
   module swap PrgEnv-intel PrgEnv-cray
   module load cray-netcdf cray-hdf5
   ./maketools -n GRIDGEN -m XC_ARCHER
   module swap PrgEnv-cray PrgEnv-intel
   export TDIR=$WDIR/dev_r4621_NOC4_BDY_VERT_INTERP/NEMOGCM/TOOLS

.. note:: my standard ARCHER ENV is intel with parallel netcdf you may need to edit accordingly

Back in $WDIR/INPUTS, create a new coordinates file from the
existing global 1/12 mesh and refine to 1/84 degree resolution: 

::
 
   cd $TDIR/GRIDGEN
   cp $WDIR/INPUTS/namelist_R12 ./
   ln -s namelist_R12 namelist.input
   ./create_coordinates.exe 
   cp 1_coordinates_ORCA_R12.nc $WDIR/INPUTS/coordinates.nc
   module swap PrgEnv-cray PrgEnv-intel
   module unload cray-netcdf cray-hdf5
   module load cray-netcdf-hdf5parallel cray-hdf5-parallel

To create the bathymetry we use the gebco dataset. On ARCHER I
had to use a non-default nco module for netcdf operations to work.
I also had to cut down the gebco data as the SCRIP routines failed
for some unknown reason.

::

   cd $WDIR/INPUTS
   module load nco/4.5.0
   ncap2 -s 'where(topo > 0) topo=0' gebco_1_cd_v2.nc tmp.nc
   ncflint --fix_rec_crd -w -1.0,0.0 tmp.nc tmp.nc gebco_in.nc
   rm tmp.nc
   $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_gebco
   $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_gebco
   $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_gebco
    
We perform a similar operation to create the initial conditions:

.. note:: I've put a sosie pre-step in here to flood fill the land. 
          I tried using sosie for 3D intepolation, but not convinced.

::

   cd ~
   mkdir local 
   svn co svn://svn.code.sf.net/p/sosie/code/trunk sosie
   cd sosie
   cp $WDIR/INPUTS/make.macro ./
   make
   make install
   export PATH=~/local/bin:$PATH   
   cd $WDIR/INPUTS
   sosie.x -f initcd_votemper.namelist
   sosie.x -f initcd_vosaline.namelist
   $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_initcd_votemper
   $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_initcd_votemper
   $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_votemper
   $TDIR/WEIGHTS/scripinterp.exe namelist_reshape_bilin_initcd_vosaline

Finally we setup weights files for the atmospheric forcing:

::

   $TDIR/WEIGHTS/scripgrid.exe namelist_reshape_bilin_atmos
   $TDIR/WEIGHTS/scrip.exe namelist_reshape_bilin_atmos
   $TDIR/WEIGHTS/scripshape.exe namelist_reshape_bilin_atmos


Next step is to create the mesh and mask files that will be used 
in the generation of the open boundary conditions:

::

   cd $CDIR
   cp $WDIR/INPUTS/cpp_LH_REEF.fcm LH_REEF/
   ln -s $WDIR/INPUTS/bathy_meter.nc $CDIR/LH_REEF/EXP00/bathy_meter.nc 
   ln -s $WDIR/INPUTS/coordinates.nc $CDIR/LH_REEF/EXP00/coordinates.nc 
   cp $WDIR/INPUTS/runscript $CDIR/LH_REEF/EXP00
   cp $WDIR/INPUTS/namelist_cfg $CDIR/LH_REEF/namelist_cfg
   cp $WDIR/INPUTS/namelist_ref $CDIR/LH_REEF/namelist_ref
   ./makenemo -n LH_REEF -m XC_ARCHER_INTEL
   cd LH_REEF/EXP00
   ln -s /work/n01/n01/jdha/ST/xios-1.0/bin/xios_server.exe xios_server.exe
   qsub -q short runscript

.. note:: there is an assumption that XIOS is already installed

If that works, we then need to rebuild the mesh and mask files in 
to single files for the next step:

::

   $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_zgr 96
   $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mesh_hgr 96
   $TDIR/REBUILD_NEMO/rebuild_nemo -t 24 mask 96
   mv mesh_zgr.nc mesh_hgr.nc mask.nc $WDIR/INPUTS
   cd !$
   rm mesh_* mask_*

Now we're ready to generate the boundary conditions using pyNEMO. 
If this is not installed follow the `installation guide` or a quick
setup could be as follows:

:: 

   cd ~
   conda create --name pynemo python scipy numpy matplotlib basemap netcdf4   
   soure activate pynemo
   conda install -c https://conda.anaconda.org/srikanthnagella seawater
   conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
   conda install -c https://conda.anaconda.org/srikanthnagella pyjnius
   export LD_LIBRARY_PATH=/opt/java/jdk1.7.0_45/jre/lib/amd64/server:$LD_LIBRARY_PATH
   export PYTHONPATH=/home/n01/n01/jdha/.conda/envs/pynemo/lib/python2.7:$PYTHONPATH
   export PYTHONPATH=/home/n01/n01/jdha/.conda/envs/pynemo/lib/python2.7/site-packages/:$PYTHONPATH
   export PYTHONPATH=/home/n01/n01/jdha/.conda/envs/pynemo/lib/python2.7/site-packages/lib/python2.7/site-packages/:$PYTHONPATH
   export PATH=/home/n01/n01/jdha/.conda/envs/pynemo/lib/python2.7/site-packages/bin/:$PATH
   svn checkout http://ccpforge.cse.rl.ac.uk/svn/pynemo
   cd pynemo/trunk/Python
   python setup.py build
   python setup.py install --prefix /home/n01/n01/jdha/.conda/envs/pynemo/lib/python2.7/site-packages/
   cd $WDIR/INPUTS

Start up pynemo and generate boundary conditions. First we need to
create a few ncml files to gather input data and map variable names.
Then using pynemo we define the area we want to model:

::

   pynemo_ncml_generator   
   pynemo -g -s namelist.bdy

Let's have a go at running the model:

::

   cd $CDIR/LH_REEF/EXP00
   sed -e 's/nn_msh      =    3/nn_msh      =    0/' namelist_cfg > tmp
   sed -e 's/nn_itend    =       1 /nn_itend    =       120 /' tmp > namelist_cfg
   qsub -q short runscript
