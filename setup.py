from distutils.core import setup
import setuptools
import os

setup(
    name='GenomeRunnerSNP',
    version='0.1.0',
    author='Mikhail Dozmorov',
    author_email='mikhail-dozmorov@omrf.org',
    packages=['grsnp'],
    package_dir={"grsnp": "grsnp"},
    scripts=[os.path.join(r,f) for r,d,fs in os.walk("grsnp/bin") for f in fs if f.endswith(".py") or f.endswith(".sh") or "bedToBigBed" in f],
    url='',
    license='LICENSE.txt',
    description='',
    long_description=open('README').read())