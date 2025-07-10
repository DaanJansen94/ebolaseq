from setuptools import setup, find_packages

setup(
    name="ebolaseq",
    version="0.1.4",
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
