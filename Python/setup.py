from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'DESCRIPTION.rst'),encoding='utf-8') as f:
    long_description = f.read()
    
setup(
      name = 'pynemo',
      
      version='0.2',
      
      description = 'NEMO Regional Configuration Toolbox',
      long_description = long_description,
      
      #The projectt's main homepage
      url = 'http://ccpforge.cse.rl.ac.uk/gf/project/src/',
      
      #Author details
      author='James Harle, Srikanth Nagella, Shirley Crompton',
      author_email='srikanth.nagella@stfc.ac.uk',
      
      #Choose your license
      license='GPL',
      
      classifiers=[
                   'Development Status :: 3 - Alpha',
                   
                   'Intended Audience :: Developers',
                   'Topic :: Oceonography Modelling',
                   'License :: OSI Approved :: GPL License',
                   
                   #Specify the python versions supported
                   'Programming Language :: Python :: 2.7'
                   ],
      
      keywords='Oceanography NEMO',
      
      packages=['pynemo','pynemo.tests','pynemo.gui','pynemo.utils','pynemo.tide','pynemo.reader'],
      
      install_requires=['netCDF4','scipy','numpy','matplotlib', 'basemap', 'thredds_crawler'],
      
      include_package_data=True,
      #The data files that needs to be included in packaging
      package_data={'': ['gui/*.png','reader/jars/*.jar'],
                    
                    },
      #If files are needs outside the installed packages
      data_files=[],
      
      entry_points={
                    'console_scripts':[
                        'pynemo=pynemo.pynemo_exe:main',
                        'pynemo_settings_editor=pynemo.pynemo_settings_editor:main'
                        ],
                    },
      )

