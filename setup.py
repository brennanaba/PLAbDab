from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='PLAbDab',
    version='0.0.1',
    description='Set of functions to search PLAbDab',
    license='BSD 3-clause license',
    maintainer='Brennan Abanades',
    long_description=long_description,
    long_description_content_type='text/markdown',
    maintainer_email='brennan.abanadeskenyon@stx.ox.ac.uk',
    include_package_data=True,
    packages=find_packages(include=('PLAbDab', 'PLAbDab.*')),
    install_requires=[
        'numpy',
        'numba',
        'kasearch',
        'pandas',
    ],
)