#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='PYMEcs',
      version='1.1',
      description='Extra functionality for PYME by CS',
      author='Christian Soeller',
      author_email='c.soeller@gmail.com',
      packages=find_packages(),
      package_data={
            # include all svg and html files, otherwise conda will miss them
            # '': ['*.svg', '*.html'],
      },
      entry_points = {
        'console_scripts': ['PYMEconfigutils=PYMEcs.misc.configUtils:main'],
      }
     )
