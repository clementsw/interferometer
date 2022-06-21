# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 16:16:20 2018

@author: clementsw
"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='interferometer',
      version='1.0',
      description='Algorithms for universal interferometers',
      long_description='Decompose unitary matrices into interferometers using two different algorithms, and visualize the interferometers',
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Physics',
      ],      
      url='https://github.com/clementsw/interferometer',
      author='William R. Clements',
      author_email='mail@william-clements.com',
      license='MIT',
      packages=['interferometer'],
      install_requires=[
          'numpy',
          'matplotlib',
      ],
      extras_require={
        'tests': ['pytest'],
        'extra': ['jupyter']
      },
      include_package_data=True,
      zip_safe=False)
