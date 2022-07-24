from setuptools import setup, find_packages

setup(
    name='elastic_vasp',
    version='0.0.1',    
    description='A example Python package',
    url='https://github.com/prnvrvs/elastic_vasp',
    author='PRANAV KUMAR',
    author_email='prnvkmr4@gmail.com',
    license='MIT Licence',
    packages=['elastic_vasp'],
    install_requires=['numpy',
                      'matplotlib'                     
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)

