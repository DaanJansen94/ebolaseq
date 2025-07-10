from setuptools import setup, find_packages

# Try to get version from setuptools_scm, fallback to hardcoded version
try:
    from setuptools_scm import get_version
    version = get_version()
except (ImportError, LookupError):
    # Fallback version when setuptools_scm fails (e.g., in conda-build)
    version = "0.1.5"

setup(
    name="ebolaseq",
    version=version,
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