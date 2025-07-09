from setuptools import setup, find_packages

setup(
    name="ebolaseq",
    use_scm_version={
        "version_scheme": "post-release",
        "local_scheme": "dirty-tag"
    },
    setup_requires=["setuptools_scm"],
    description="A tool for filtering and analyzing Ebola virus sequences",
    author="Daan Jansen",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "biopython>=1.81",
        "setuptools_scm",
    ],
    entry_points={
        'console_scripts': [
            'ebolaseq=ebolaseq.ebolaseq:cli_main',
        ],
    },
)
