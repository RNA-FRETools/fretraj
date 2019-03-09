import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="acv_pkg",
    version="1.0.0",
    author="Fabio Steffen",
    author_email="fabio.steffen@chem.uzh.ch",
    description="Calculating accessible-contact volumes of fluorescent dyes on biomolecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fdsteffen/acv",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords='ACV, FRET, PDB'
)
