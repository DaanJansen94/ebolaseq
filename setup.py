from setuptools import setup, find_packages

# Import version from the package
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'ebolaseq'))
from __init__ import __version__

setup(
    name="ebolaseq",
    version=__version__,
    description="A tool for filtering and analyzing Ebola virus sequences",
    author="Daan Jansen",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "biopython>=1.81",
    ],
    entry_points={
        'console_scripts': [
            'ebolaseq=ebolaseq.ebolaseq:cli_main',
        ],
    },
)
