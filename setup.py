from setuptools import setup, find_packages

requires = [
    'biopython>=1.7',
    'pandas'
    ]

setup(
    name="bioutils",
    version="1.0",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("Bioinformatic toolkit for processing and analysing omic data"),
    license="LICENSE",
    #install_requires=requires,
    packages=['bioutils'],
    entry_points={
        'console_scripts': ['bioutils=bioutils.main:main'],
    }
)
