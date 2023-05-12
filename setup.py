from setuptools import setup, find_packages

VERSION = '1.0'
DESCRIPTION = 'BGCs screening for antibacterial resistence hits'
LONG_DESCRIPTION = 'Biosynthetic gene clusters (BGCs) screening for antibacterial resistence hits, aiming to prospect new possible antibiotics'

setup(
	name='Krill',
	version=VERSION,
	author='Saulo Britto da Silva',
	author_email='<saulobdasilva@gmail.com>',
	description=DESCRIPTION,
	long_description=LONG_DESCRIPTION,
	url='',
	license='',
	packages=find_packages(),
	install_requires=['tqdm','cprint']
	)

keywords=['BGC','screening','gene clusters','biosynthetic','secondary metabolities']
classifiers=[
			'Development Status :: 2 - Pre-Alpha',
			'Intended Audience :: Science/Research',
			'License :: ',
			'Programming Language :: Python :: 3 :: Only',
			'Operating System :: POSIX :: Linux'
			]