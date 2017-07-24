#SAAV-structural-mapping
![](http://i.imgur.com/NSh2bGe.gif)

### Description

The contents of this repository are being continually developed to understand how genetic variation within a population influences the structure and function of proteins, and importantly, vice versa.

### Requirements

1. PyMOL
	* All molecular visualizations are done using PyMOL (https://pymol.org). If you do not have a licensed version, not all features will be available to you. With the educational version everything will work except that the output images will be low quality (https://pymol.org/edu/?q=educational). With the open source version (https://sourceforge.net/projects/pymol/) SAAV coloring, transparency, and radius size will not work as expected :(


2. Python 2
	* All python programs with the line `import pymol` are written for Python 2. This is because it makes use of the program PyMOL (https://pymol.org/) for molecular visualization which is only developed for Python 2. If you want to go down the rabbit hole, you should be able to run the open source PyMOL version on Python 3 (https://sourceforge.net/projects/pymol/), but any of the programs with `import pymol` will need to be converted. Also, unfortunately some features we utilize are not available in the open source code :(

### Usage

