# *f*-distance of knotoids

This program computes the experimental f-distance of isotopy classes of knotoid diagrams, as presented in [1]. We use the classification data for knotoids on the sphere from [2].


## Dependencies
The script has been tested on Python 3.7 on OSX 10.14.5 and on Ubuntu Linux 16.04.6. The following python packages are required:

1. pandas
2. python-igraph

To see the installed packages of your python distribution type <code>pip list</code>. If you use conda type <code>conda list</code>.
To install pandas using Pip type <code> pip install pandas </code>. If you use conda type <code> conda install -c anaconda pandas </code>.

To install python-igraph on OSX use the following commands in the terminal:
<code> brew install cairo </code>,
<code> brew install py2cairo </code>,
<code> brew install igraph  </code> and
<code> pip install python-igraph </code> for pip users or <code>conda install -c conda-forge python-igraph</code> for conda users.

Note that for the above commands to work, Homebrew is required. For installation please follow the instructions at: http://brew.sh/ 
Alternatively on may use Macports: https://www.macports.org/install.php

To install python-igraph on Linux, use the command <code> apt-get install python-igraph </code>



## Usage

To use the program one has first to load a file containing knotoid data as follows: `python f-distance.py -f INPUTFILENAME`. The output file will have the name: `INPUTFILENAME_forbidden_data.txt`

The program computes by default the f-distance only for *prime* knotoids. To add composites into the analysis use the option -c: `python bracket.py -f INPUTFILENAME -c`. 

To plot the resulting graph use the option -p: `python f-distance.py -f INPUTFILENAME -c -p`. 


## Input file format

The input file must be a tab separated file with the following columns:
gc1,	gc2,	(optional: gc3),	status,	isotopy.class

gc1: list of crossings labels as one walks around the diagram. An overcrossing is noted with a positive sign while an undercrossing with a negative sign.

gc2: list of crossing signs. Right handed crossing is positive, left handed is negative.

gc3: OPTIONAL. Used only in the case of planar knotoids.

status: status of the diagram. Can be prime, not_prime or composite.

isotopy.class: The isotopy class of the diagram.

The file `primes_mirror_symmetric_sphere.txt` is included as an example. It is taken from [2].

## Credits

### Contributors
Dimos Goundaroulis

### Institutions
Center for Integrative Genomics, University of Lausanne.

SIB Swiss Institute of Bioinformatics.


## References
[1] A. Barbensi and D. Goundaroulis *f-distance of knotoids and protein structure*. arXiv preprint (2019).
 
[2] D. Goundaroulis, J. Dorier and A. Stasiak, *A systematic classification of knotoids on the plane and on the sphere*. arXiv preprint [arXiv:1902.07277](https://arxiv.org/abs/1902.07277) (2019).
