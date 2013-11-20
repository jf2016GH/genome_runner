from distutils.core import setup
import setuptools
import os

# Create list of data files with paths relative to the base genome-runner directory
package_data = []
for root, dirs, files in os.walk("grsnp/frontend"):
    for file in files:
        package_data.append(os.path.relpath(os.path.join(root, file), "grsnp/frontend"))

setup(
    name='GenomeRunner SNP',
    version='0.1.0',
    author='Mikhail G. Dozmorov, Lukas R. Cara, Cory B. Giles',
    author_email='mikhail.dozmorov@gmail.com, lks_cara@yahoo.com, mail@corygil.es',
    data_files=[('grsnp', ['grsnp_db_readme.txt']),],
    packages=['grsnp'],
    package_dir={"grsnp": "grsnp"},
    package_data={"": package_data},
    scripts=[os.path.join(r,f) for r,d,fs in os.walk("grsnp/bin") for f in fs if f.endswith(".py") or f.endswith(".sh") or "bedToBigBed" in f],
    url='http://www.genomerunner.org',
    license='LICENSE.txt',
    description='GenomeRunner SNP: Interpreting genome veriation within epigenomic context',
    long_description=open('README').read())
