"""tool"""
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


install_requires = ['setuptools>=18.0', 'MOODS-python']


setup(
    name='cenrom',
    version='0.0.2',
    description='calculate enrichment of matrix in data',
    author='Anton Tsukanov',
    author_email='tsukanov@bionet.nsc.ru',
    url='http://github.com/ubercomrade/pipeline',
    package_dir={'lib' : 'lib'},
    packages=[
        'lib',
    ],
    scripts=['cenrom.py',],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        'Programming Language :: Cython',
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe=False,
    include_package_data=True,
    install_requires=install_requires,
    setup_requires=install_requires,
    python_requires='>=3.7',
)