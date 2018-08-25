![alt text](https://github.com/btozer/gmg/blob/master/gmg/docs/_sources/_static/gmg_logo.png)


|PyPI| |Affiliated package| |Coverage Status| |Build status|

Statement of need
=================

GMG is an open-source Graphical User Interface (GUI) designed principally for modelling
2D potential field (gravity and magnetic) profiles. The software also includes 
functions for loading XY data, seismic reflection SEGY data and exploration well horizons.
The software therefore provides an integrated geological/geophysical interpretation
package. It is anticipated that GMG will also be useful for teaching purposes.

Data I/O is made as simple as possible using space delimited ASCII text files.

The project was instigated after failing to find an adequate open source option
(in which the source code can be viewed and modified by the user) for performing 2D 
geophysical modeling tasks. Inspiration came from fatiando a terra and GMT.

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
* Display XY data (e.g., earthquake hypocenters or geological surface contacts).
* Export model data (e.g., predicted anomalies and layer geometries) as ASCII text files.
* Save model figures as vector or raster graphics in various formats (.png .ps .eps .pdf).


Installation
============

Pre-Installation
----------------

To install gmg you will need some familiarity with using a unix terminal. The simplest way to install gmg is to use an Anaconda Python distribution. Anaconda can be downloaded from [https://www.continuum.io/downloads](https://www.continuum.io/downloads)

Using Anaconda will ensure you can install all the dependencies required on any platform (Linux, Mac, Windows) 
using the Conda package manager. Additionally, you can create a separate conda environment to use when running gmg. 
This will avoid any potential conflicts with your current system configuration.

To create a new python2.7 conda environment you can use the command:

    conda create -n py27-gmg python=2.7 anaconda

Where the -n flag dictates what you want to name the new environment. 
This environment can then be activated using:

    source activate py27-gmg

Supported OS
------------

**NB:** gmg is (for the time being) only tested on **Ubuntu 14.04 Linux** and **macOS High Sierra** operating systems. 
It is anticipated that issues may arise when trying to install on other distributions, particularly Windows.
Please raise any information regarding installation problems on github.


Installing gmg
--------------

**Step 1:**

Ensure your py27-gmg environment is active:

    source activate py27-gmg


**Step 2: Install dependencies**

GMG depends on several other packages to run. These are:

* [wxpython](http://wiki.wxpython.org/)
* [numpy](http://www.numpy.org)
* [scipy](http://scipy.org/)
* [matplotlib](http://matplotlib.sourceforge.net/)
* [fatiando a terra](http://www.fatiando.org/)
* [ObsPy](http://docs.obspy.org/)

These dependencies can be installed in your py27 environment using the command:

    conda install wxpython numpy scipy matplotlib fatiando obspy


**Step 3: Get gmg**

Download or git clone the gmg github repository. To download, navigate to:

    https://github.com/btozer/gmg

Use the green *"clone or download"* button on the right hand side of the page to download a .zip of gmg. Unzip this
directory.

or simply use git on the command line:

    git clone https://github.com/btozer/gmg.git


**Step 4: Install  gmg**

Run the following command in the dir where you have downloaded or git-cloned the gmg repository:

    pip install .

This will install gmg within your local anaconda python site-packages dir.


Launching gmg
-------------

To start gmg the best way is to create an alias and export this, e.g.:

    export alias gmg='python ~/anaconda3/envs/py27-gmg/lib/python2.7/site-packages/gmg/gmg.py'

or manually add this to your .bashrc or equivalent file (e.g. .bash_profile, .zshrc or .cshrc) e.g.:

    alias gmg='python ~/anaconda3/envs/py27-gmg/lib/python2.7/site-packages/gmg/gmg.py'

NB. for macOS users you may need to invoke pythonw instead of python if you receive the message:

    This program needs access to the screen. Please run with a
    Framework build of python, and only when you are logged in
    on the main display of your Mac.

i.e.:

    alias gmg='pythonw ~/anaconda3/envs/py27-gmg/lib/python2.7/site-packages/gmg/gmg.py'

Then you can (optionally) add the executable directory to your $PATH variable by repeating the process, e.g.:

    export $PATH=$PATH:~/anaconda3/envs/py27-gmg/lib/python2.7/site-packages/gmg/

Now, in a fresh terminal, simple type "gmg" on the command line to launch the software:

    gmg


Getting started
---------------

It is recommended you run the two test cases when first launching the software to check the potential field algorithms 
are running correctly. Details can be found in the documentation. This can be accessed from within
the software under:

    Help -> Documentation

or by opening:
    
    "PATH_TO_GMG"/gmg/docs/_build/html/gmg_documentation.html

If the test cases don't work as expected, please submit a bug report as described
below in the section *Contributing*.

Once the test cases have been run (and work as expected) the easiest way to get a feel for the 
software is to work through the tutorial (see *documentation*).


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
