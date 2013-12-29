from distutils.core import setup
from distutils.command.install import install

import setuptools
import os
import suprocess

# Create list of data files with paths relative to the base genome-runner directory
package_data = ["frontend/*"]

def scrip_installer(command_subclass):
    """A decorator for classes subclassing one of the setuptools commands.

    It modifies the run() method so that is runs scripts to install required packages.
    Only works in Ubuntu
    """
    orig_run = command_subclass.run

   def modified_run(self):
        # Install the R packages required by grsnp
        r_packages_install = "wget http://cran.r-project.org/src/contrib/Hmisc_3.13-0.tar.gz\nwget http://cran.r-project.org/src/contrib/gplots_2.12.1.tar.gz\nwget http://cran.r-project.org/src/contrib/RColorBrewer_1.0-5.tar.gz\nsudo R CMD INSTALL Hmisc_3.13-0.tar.gz gplots_2.12.1.tar.gz RColorBrewer_1.0-5.tar.gz\nrm Hmisc_3.13-0.tar.gz gplots_2.12.1.tar.gz RColorBrewer_1.0-5.tar.gz"
        out = subprocess.Popen(r_packages_install,stdout=subprocess.PIP,shell=True)
        out.wait()
        # installs grtk
        grtk_install = """curl http://bedops.googlecode.com/files/bedops_linux_x86_64-v2.2.0.tar.bz2 | sudo tar xvj -C /usr/local\nsudo wget -np -R -A "bedToBigBed" -A "bedGraphToBigWig" -A "bigWig*" -A "bigBed*" -N -e robots=off -r -P /usr/local/bin -nd "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"\nsudo curl -o /usr/local/bin/rowsToCols http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/rowsToCols\nsudo chmod a+x /usr/local/bin/*\nsudo apt-get install -y parallel bedtools tabix kyotocabinet-utils realpath\ngit clone git@bitbucket.org:wrenlab/grtk.git\ncd grtk\nsudo python setup.py install"""
        out = subprocess.Popen(grtk_install,stdout=subprocess.PIP,shell=True)
        out.wait()

        orig_run(self)

    command_subclass.run = modified_run
    return command_subclass

@friendly
class CustomInstallCommand(install):
    pass

setup(
    name='GenomeRunner SNP',
    version='0.1.0',
    author='Mikhail G. Dozmorov, Lukas R. Cara, Cory B. Giles',
    author_email='mikhail.dozmorov@gmail.com, lks_cara@yahoo.com, mail@corygil.es',
    data_files=[('grsnp', ['grsnp_db_readme.txt']),],
    packages=['grsnp'],
    package_dir={"grsnp": "grsnp"},
    package_data={"grsnp": package_data},
    scripts=[os.path.join(r,f) for r,d,fs in os.walk("grsnp/bin") for f in fs if f.endswith(".py") or f.endswith(".sh") or "bedToBigBed" in f],
    url='http://www.genomerunner.org',
    license='LICENSE.txt',
    install_requires=["cherrypy","numpy","scipy","cython","pybedtools","bx-python","rpy2","mako","simplejson"],
    description='GenomeRunner SNP: Interpreting genome veriation within epigenomic context',
    long_description=open('README').read(),
    cmdclass={
        'install': CustomInstallCommand
    })
