from setuptools import setup
from setuptools import find_namespace_packages

# Load the README file.
with open(file="README.md", mode="r") as readme_handle:
    long_description = readme_handle.read()

setup(
    name='AutoConverge',
    author='Kian Pu',
    author_email='jiankunp@andrew.cmu.edu',
    # Define the version of this library.
    # Read this as
    #   - MAJOR VERSION 0
    #   - MINOR VERSION 1
    #   - MAINTENANCE VERSION 0
    version='0.1.0',
    description='A python library to automate convergence test in GPAW',
    url='https://github.com/kianpu34593/AutoConverge',
    license='MIT',
    install_requires=[
        "numpy==1.19.5",
        "ase==3.20.1",
        "gpaw==20.10.0",
        "autocat==0.0.1",
        "pymatgen==2020.12.31",
        "pubchempy==1.0.4",
    ],
    keywords='DFT,GPAW,convergence test,ase',
    #package_dir={"": "GPAW_converge"},
    packages=find_namespace_packages(exclude='example'),
    python_requies='>=3.7',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.8",
    ],

)