import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

about = {}
with open('fretraj/__about__.py') as a:
    exec(a.read(), about)

setuptools.setup(
    name=about['__title__'],
    version=about['__version__'],
    author=about['__author__'],
    author_email=about['__email__'],
    description=about['__description__'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    url=about['__url__'],
    packages=setuptools.find_packages(exclude=['docs', 'tests']),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    keywords=about['__keywords__']
)
