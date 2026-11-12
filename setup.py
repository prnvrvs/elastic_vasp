import os
import re

from setuptools import setup

# Single source of truth for the version: elastic_vasp/__init__.py
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'elastic_vasp', '__init__.py')) as fh:
    VERSION = re.search(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', fh.read(), re.M).group(1)

setup(
    name='elastic_vasp',
    version=VERSION,
    description='Second-order elastic constants from VASP energy-strain calculations',
    url='https://github.com/prnvrvs/elastic_vasp',
    author='PRANAV KUMAR',
    author_email='prnvkmr4@gmail.com',
    license='MIT Licence',
    packages=['elastic_vasp'],
    python_requires='>=3.9',
    install_requires=['numpy',
                      'matplotlib',
                      # POSCAR reading/writing; older ASE ignores negative and
                      # per-axis (VASP 6) scale factors
                      'ase>=3.23',
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)

