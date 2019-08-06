![alt text](https://github.com/btozer/gmg/blob/master/docs/_sources/_static/gmg_logo.png)

---
[Visit the GMG documentation page](https://btozer.github.io/gmg/)
---

|PyPI| |Affiliated package| |Coverage Status| |Build status|


Statement of need
=================

GMG is an open-source Graphical User Interface (GUI) designed principally for modelling
2D potential field (gravity and magnetic) profiles. The software also includes 
functions for loading XY data, seismic reflection SEGY data and exploration well horizons.
The software therefore provides an integrated geological/geophysical interpretation
package. It is anticipated that GMG will also be useful for teaching purposes.

Data I/O is made as simple as possible using space delimited ASCII text files.

The project was instigated after failing to find an adequate open-source option
(in which the source code can be viewed and modified by the user) for performing 2D 
geophysical modeling tasks. Inspiration came from fatiando a terra and GMT.

**NB: GMG is in development. Some documentation is incomplete and some features may not work as expected.**

Key features
------------

* Import and display observed topography, gravity and magnetic 2D profile data.
* Apply filters to observed data.
* Import and display seismic reflection data.
* Add and manipulate model layers (subsurface bodies) using a simple interactive interface.
* “Pinch out/snap” layers against adjacent layers.
* Calculate the predicted gravity anomaly produced by any combination of model layers.
* Calculate the predicted magnetic anomaly produced by any combination of model layers.
* Model magnetic anomalies using induced and/or remanent magnetism.
* Display well horizon tops.
* Display XY data (e.g., earthquake hypocenters, geological surface contacts or rays from a auxiliary velocity model).
* Export model data (e.g., predicted anomalies and layer geometries) as ASCII text files.
* Save model figures as vector or raster graphics in various formats (.png .ps .eps .pdf).


Installation
============

Pre-Installation
----------------

**Step 1: Install an Anaconda python distribution**

The simplest way to install GMG is to first install an Anaconda Python distribution: www.anaconda.com/download

This should ensure you can install all the dependencies required on any platform (Linux, Mac, Windows)
using the Conda package manager.  GMG is written in python3 so you will need the python3 version of Anaconda.

**Step 2: Create a new python environment for gmg**

It may be useful to create a new conda environment to use when running gmg. This will avoid any potential
conflicts with your current system configuration. e.g.:

    conda create -n gmg-env python=3.7 anaconda

Where the -n flag dictates what you want to name the new environment (you can call the environment whatever you like,
but the documentation will use the name gmg-env).

On macOS and Linux this environment can then be activated using:

    source activate gmg-env

On Windows use:

    activate gmg-env

Once the environment is activated any call to python will only "see" the packages installed within the py27-gmg
environment. To deactivate the environment either close the terminal window or on macOS and Linux run:

    source deactivate

On Windows use:

    deactivate

**OS Support:**

GMG is (for the time being) only tested on **Ubuntu 14.04 Linux** and **macOS Mojave** operating systems. 
It is possible that issues may arise when trying to install on other distributions. Please raise any information 
regarding installation problems on github.



Installing gmg
--------------

**Step 1:**

Ensure your gmg-env environment is active, on macOS or Linux use:

    source activate gmg-env

Or on Windows use:

    activate gmg-env

**Step 2: Install dependencies**

GMG depends on several other packages to run. These are:

* [wxpython](http://wiki.wxpython.org/)
* [numpy](http://www.numpy.org)
* [scipy](http://scipy.org/)
* [matplotlib](http://matplotlib.sourceforge.net/)
* [fatiando a terra](http://www.fatiando.org/)
* [ObsPy](http://docs.obspy.org/)

These dependencies can be installed in your gmg-env environment using the command:

    conda config --add channels conda-forge
    conda install wxpython numpy scipy matplotlib fatiando obspy

**Step 3: Get gmg**

Download or "git clone" the gmg github repository. 

**Option 1. Download:** 

In your browser, navigate to:

    https://github.com/btozer/gmg

Use the green *"clone or download"* button on the right hand side of the page to download a .zip of GMG. Unzip this
directory.

**Option 2. git clone:**

If you have git installed, simply use git on the command line:

    git clone https://github.com/btozer/gmg.git

This will create a copy of the GMG repository in the current working directory. 

**Step 4: Install  gmg**

Run the following command in the gmg/ root directory (which contains the script **setup.py**):

    pip install .

This will install gmg within your local anaconda python site-packages dir.


Launching gmg
-------------

To start gmg, the best way is to create an alias within your .bashrc or equivalent file 
(e.g. .bash_profile, .bashrc, .zshrc or .cshrc) e.g.:

    vi .bashrc

and then add the line:

    alias gmg='python ~/anaconda3/envs/py27-gmg/lib/python2.7/site-packages/gmg/gmg.py'

NB. for macOS users you may need to invoke pythonw instead of python if you receive the message:

    This program needs access to the screen. Please run with a
    Framework build of python, and only when you are logged in
    on the main display of your Mac.

i.e.:

    alias gmg='pythonw ~/anaconda3/envs/py27-gmg/lib/python2.7/site-packages/gmg/gmg.py'

Now, in a fresh terminal, simple type "gmg" on the command line to launch the software:

    gmg


Getting started
---------------

The best way to get started is to read through the Manual and then work your way through the Tutorial. 
These can be found by launching gmg and then navigating to:

    Help -> Documentation

or by opening:
    
    "PATH_TO_GMG"/gmg/docs/_build/html/gmg_documentation.html

Contributing
------------

Any contribution, big or small is welcome. Examples include:

1. Writing new code/functions.
2. Improving the user interface.
3. Submitting bug reports and spelling corrections.

They are all helpful.

**Becoming a developer**

The best way to become a developer is to create a folk of the github repository.
Please see the documentation for details.

**Suggestions & bug reporting**

The best way to open a bug report is to raise an issue on Github.

If you don't use Github, you can send a bug report to: btozer@ucsd.edu

Please use the subject line:
    
    GMG bug report: "short bug description"

and include as much detail as possible when describing the bug.

Citing:
-------

If you use GMG, please cite this github repository.
