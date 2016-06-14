#!/usr/bin/env python

from distutils.core import setup

setup(name='prosic',
      version='1.0',
      description='A caller for somatic insertions and deletions.',
      author='Louis Dijkstra',
      url='https://github.com/prosic/prosic',
      packages=['prosic'],
      scripts=['bin/prosic-extract-observations', 'bin/prosic-annotate'])
