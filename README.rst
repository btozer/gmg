GMG: An open source two-dimensional geophysical modelling GUI
=============================================================

|PyPI| |Affiliated package| |Coverage Status| |Build status|

Statement of need
-----------------

GMG is an Open Source Graphical User Interface (GUI) designed principally for modelling
2D potential field (gravity and magnetic) profiles. The software also includes 
functions for loading XY data, seismic reflection SEGY data and exploration well horizons.
The software therefore provides an integrated geological/geophysical interpretation
package. It is anticipated that GMG will also be useful for teaching purposes.

Input and output of data is made as simple as possible using space delimited
ASCII text files.

The project was instigated after failing to find an adequate open source option
(in which the source code can be viewed and modified by the user) for performing 2D 
geophysical modeling tasks. Inspiration came from fatiando a terra and GMT.


Installation and Dependencies
------------------------------
GMG is written in Python 2.7 and relies on the dependencies listed below. 
All of these dependencies need to be installed in order for GMG to run.


**Dependencies**

    python2.7
    fatiando
    matplotlib
    numpy
    obspy
    scipy
    wxpython


**Installing**

The simplest way to install gmg is to first install an Anaconda Python 
distribution: https://www.continuum.io/downloads
This should ensure you can install all the dependencies required  on 
any platform (Linux, Mac, Windows) using the Conda package manager.

Then gmg may be installed using pip (https://packaging.python.org/tutorials/installing-packages/).

The following command will install GMG and all of its dependencies 
within your local anaconda python site-packages dir.

Using pip all of the gmg dependencies can be installed.

This will install and all of its dependencies with your local
anaconda python site-packages directory.

gmg can then be installed by running the following command in the dir where you have the source code downloaded. 
This will install gmg into your local anaconda python site-packages dir.

pip install .

or for a local installation

pip install .

**Launching gmg**

To start gmg the best way is to create an alias and export this e.g.:

export alias gmg='python ~/anaconda/lib/python3/site-packages/gmg/gmg.py'

or manually add this to your .bashrc or equivalent file (e.g. .bash_profile or .cshrc) e.g.:
alias gmg='python ~/anaconda/lib/python3/site-packages/gmg/gmg.py'

Then add the executable directory to your $PATH variable by repeating the progess, e.g.:
export $PATH=$PATH:~/anaconda/lib/python3/site-packages/gmg/

Now, in a fresh terminal, simple type "gmg" on the command line to launch the software.


Getting started
---------------

It is recommended you run the two test cases when first launching the software
to check the potential field algorithms are running correctly.
Details can be found in the documentation. This can be accessed from within
the software under:

Help -> Documentation

or

by opening PATH_TO_GMG/gmg/gmg/docs/Manual.html

If the test cases don't work as expected, please submit a bug report as described
below in the section **Contributing**.

Once the test cases have been run (and work) the easiest way to get a feel for the 
software is run through the tutorial as described in the documentation.


Contributing
------------

Any contribution, big or small is welcome. Examples include:

1. Writing new code/functions.
2. Improving the user interface.
3. Making suggestions for improvements or brand new functions.
4. Submitting bug reports and spelling corrections.

They are all helpful.

**Becoming a developer**

The best way to become a developer is to create a folk of the github repository.
Please see the documentation for details.

**Suggestions & bug reporting**

The best way to open a bug report is to raise an issue on Github.

If you don't use Github, you can send a bug report to: btozer@ucsd.edu
Please use the subject line - GMG bug report: "short bug description"
and include as much detail as possible when describing the bug.
