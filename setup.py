# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:08:37 2019

@author: Colin Helms
"""

from distutils.core import setup
setup(name='alfano',
      version='0.4a',
      long_description='Libraries for trajectory solutions to Edelbaum-Alfano low-thrust orbit transfer.',
      description='Edelbaum-Alfano libraries.',
      author='Colin Helms',
      author_email='colinhelms@outlook.com',
      url='https://www.FreelanceRocketScience.com',
      packages=['', 'alfano', 'alfano.utilities', 'alfano.controls'],
      package_data={'alfano.controls' : ['Controls.json']},
      license='EULA',
      classifiers=[
              'Development Status :: 2 - Alpha',
              'Environment :: Console',
              'Intended Audience :: End Users',
              'Operating System :: Microsoft :: Windows7',
              'Operating System :: Linux :: POSIX'
              'Programming Language :: Python :: 3.5',
              'License :: End User License Agreement',
              'Topic :: Copyright :: Copyright Freelance Rocket Science, 2019']
      )