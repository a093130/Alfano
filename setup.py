# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:08:37 2019

@author: Colin
"""

from distutils.core import setup
setup(name='alfano',
      version='0.2a',
      author='Colin Helms',
      author_email='colinhelms@outlook.com',
      url='https://www.FreelanceRocketScience.com/downloads',
      packages=['utilities', 'controls'],
      package_data={'controls' : ['Controls.json']},
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