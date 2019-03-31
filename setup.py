from setuptools import setup

setup(
    name='macromaker',
    version='0.1.0',
    description='A macrocomplex maker from pdb pieces',
    url='http://pypi.python.org/pypi/macromaker/',
    author='The Bioinfo Crew',
    keywords='protein pdb macrocomplex',
    packages=['macromaker', 'macromaker.lib','macromaker.test'],
    scripts=['bin/macromaker','bin/correct-pdb'],
    license='GNU GPLv3',
    long_description=open('README.md').read(),
    install_requires=[
        "biopython >= 1.73",
        "python-Levenshtein >= 0.12.0",
        "numpy >= 1.16.2",
        "pytest",
    ],
    python_requires='>=3',
)