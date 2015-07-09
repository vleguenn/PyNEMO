Installation
============
This page is installation to pyNEMO

Dependencies
^^^^^^^^^^^^

1. Python 2.7 (Not tested with 3.x)
2. scipy
3. netCDF4-python
4. numpy
5. matplotlib
6. basemap
7. thredds_crawler
8. pyjnius (optional)

Anaconda
^^^^^^^^

Using conda. pyNEMO supports Win64, OSX and Linux. for other operating systems please build from source.

::

   conda install -c https://conda.anaconda.org/srikanthnagella pynemo

This will install the pynemo and its dependencies. 

From Source
^^^^^^^^^^^

Installing pyNEMO using other flavours of software or from source. Install all the dependencies and download the source code from svn and install. 

::

   svn checkout http://ccpforge.cse.rl.ac.uk/svn/pynemo/trunk/Python/
   python setup.py install
   