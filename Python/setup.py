from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'DESCRIPTION.rst'),encoding='utf-8') as f:
    long_description = f.read()
    
setup(
      name = 'pyNEMO',
      
      version='0.1',
      
      description = 'A Python NEMO Project',
      long_description = long_description,
      
      #The projectt's main homepage
      url = 'http://ccpforge.cse.rl.ac.uk/gf/project/pynemo/',
      
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
      
      packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
      
      install_requires=['netCDF4','scipy','numpy'],
      
      #The data files that needs to be included in packaging
      package_data={
                    },
      #If files are needs outside the installed packages
      data_files=[],
      
      entry_points={
                    'console_scripts':[
                        'pyNEMO=src:profile',
                        ],
                    },
      )