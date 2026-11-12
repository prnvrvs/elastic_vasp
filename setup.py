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
    license='MIT',
    long_description=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README.md'), encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    project_urls={
        'Documentation': 'https://prnvrvs.github.io/elastic_vasp/',
        'Source': 'https://github.com/prnvrvs/elastic_vasp',
    },
    packages=['elastic_vasp'],
    python_requires='>=3.10',
    install_requires=['numpy',
                      'matplotlib',
                      # POSCAR reading/writing; older ASE ignores negative and
                      # per-axis (VASP 6) scale factors
                      'ase>=3.23',
                      # crystal symmetry for the orientation checks
                      'spglib>=2.4',
                      ],
    entry_points={
        'console_scripts': ['elastic-vasp = elastic_vasp.cli:main'],
    },
    extras_require={
        'test': ['pytest', 'scipy'],
        'doc': ['sphinx', 'sphinx-book-theme'],
    },

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
)

