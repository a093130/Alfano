#!python
# # -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:08:37 2019

@author: Colin Helms
"""

from distutils.core import setup
setup(name='alfano',
      version='0.5a1',
      description='Edelbaum-Alfano Low-Thrust Trajectories.',
      long_description='file: README.txt',
      author='Colin Helms',
      author_email='colinhelms@outlook.com',
      url='https://github.com/a093130/Alfano',
      packages=['', 'src', 'src.utilities', 'src.controls'],
      classifiers=[
              'Development Status :: Alpha',
              'Environment :: Console',
              'Operating System :: Microsoft :: Windows 10',
              'Operating System :: Linux :: Ubuntu 18.4',
              'Programming Language :: Python :: 3.4',
              'License :: End User License Agreement',
              'Topic :: Copyright :: Copyright Freelance Rocket Science, 2022'
              ],
        )