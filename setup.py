import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fretraj",
    version="1.0.0",
    author="Fabio Steffen",
    author_email="fabio.steffen@chem.uzh.ch",
    description="Calculating accessible-contact volumes of fluorescent dyes on biomolecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fdsteffen/fretraj",
    packages=setuptools.find_packages(exclude=['docs', 'tests']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    keywords='accessible volume, contact volume, MD, trajectory, ACV, FRET, PDB'
)
