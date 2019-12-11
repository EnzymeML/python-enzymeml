from setuptools import setup

setup(
    name='python-enzymeml',
    version='1.0.0',
    packages=['enzymeml'],
    url='https://github.com/EnzymeML/python-enzymeml',
    license='BSD-2-Clause',
    author='EnzymeML Team',
    author_email='enzymeml-team@googlegroups.com',
    description='Library for reading / writing Enzyme ML documents',
    install_requires=['python-libsbml', 'python-libcombine', 'xlrd'],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/EnzymeML/python-enzymeml/issues',
        'Source': 'https://github.com/EnzymeML/python-enzymeml',
    },
)
