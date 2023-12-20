from setuptools import setup, find_packages

requires = [
    'biopython>=1.7',
    'pandas'
    ]

setup(
    name="bioinf_utils",
    version="1.0",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("Bioinformatic toolkit for processing and analysing omic data"),
    license="LICENSE",
    #install_requires=requires,
    packages=['bioinf_utils'],
    entry_points={
        'console_scripts': ['bioutils=bioinf_utils.main:main'],
    }
)
