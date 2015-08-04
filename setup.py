#!/usr/bin/env python
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup

setup(name='GenomeRunner web',
	version='0.1.0',
	author='Mikhail G. Dozmorov, Lukas R. Cara, Cory B. Giles',
	author_email='mikhail.dozmorov@gmail.com, lks_cara@yahoo.com, mail@corygil.es',
	license='License.txt',
	url='http://www.genomerunner.org',
	description='GenomeRunner web: Interpretation of SNP sets within regulatory context',
	long_description=open('README.md').read(),
	packages=['grsnp'],
	package_dir={'grsnp': 'grsnp'},
	include_package_data=True, # Install data from MANIFEST.in,
	entry_points={
		'console_scripts':[
			'gr = grsnp.hypergeom4:main',
			'gr-server = grsnp.server:main',
			'gr-optimizor = grsnp.optimizor:main',
			'gr-dbcreator = grsnp.dbcreator_encode:main'
		]
	},
	install_requires = [
		"cherrypy==3.5",
		"celerly==3.1",
		"numpy==1.8",
		"scipy==0.14",
		"rpy2==2.5",
		"simplejson==3.6",
		"mako==1.0",
		"BeautifulSoup==3.2",
		"redis==2.10"
	]
)

