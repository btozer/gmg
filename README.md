![alt text](https://github.com/btozer/gmg/blob/add_vgg/docs/_sources/_static/gmg_logo_white.png)

---
[Visit the GMG documentation page](https://btozer.github.io/gmg/)
---
<p align="center">
<a href="https://test.pypi.org/project/gmgpy/">
<img
src="http://img.shields.io/pypi/v/gmgpy.svg?style=flat-square"
alt="Latest version on PyPI"
/>
</a>
<a href="https://github.com/conda-forge/gmgpy-feedstock">
<img
src="https://img.shields.io/conda/vn/conda-forge/gmgpy.svg?style=flat-square"
alt="Latest version on conda-forge"
/>
</a>


![alt-text](https://github.com/btozer/gmg/blob/add_vgg/gmg_demo.gif)

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
geophysical modeling tasks. Inspiration came from [Fatiando a Terra](https://www.fatiando.org/) and [GMT](https://www.generic-mapping-tools.org/).

**NB: GMG is in development. Some documentation is incomplete and some features may not work as expected. If you experience any issues, please feel free to raise these on the issue tracker.**

Key features
------------

* Import and display observed topography, gravity and magnetic 2D profile data.
* Apply filters to observed data.
* Add and manipulate model layers (subsurface bodies) using a simple interactive interface.
* “Pinch out/snap” layers against adjacent layers.
* Calculate the predicted gravity anomaly produced by any combination of model layers.
* Calculate the predicted magnetic anomaly produced by any combination of model layers.
* Model magnetic anomalies using induced and/or remanent magnetism.
* Display well horizon tops.
* Display XY data (e.g., earthquake hypocenters, geological surface contacts or rays from a auxiliary velocity model).
* Import and display seismic reflection data.
* Export model data (e.g., predicted anomalies and layer geometries) as ASCII text files.
* Save model figures as vector or raster graphics in various formats (.png .ps .eps .pdf).

Installation
============

Pre-Installation
----------------

**Note** 

gmg is written in Python 3 and can be installed using common python package managers 
conda (from conda-forge) or pip. Installation and launching of the software involves 
the use of a command line terminal, but this is minimal and should be manageable for a non-expert
by following the instructions below. 

**Step 1: Install Miniforge**

If you don't already use miniforge or anaconda for managing and installing python packages, then 
the simplest way to install gmg is to first install Miniforge (a minimal installer for [Conda](http://docs.conda.io/en/latest/)
specifically for [conda-forge](http://conda-forge.org)).

This should ensure you can install all the dependencies required on any platform (Linux, Mac, Windows) using the Conda
package manager.

**Step 2: Create a new conda environment for gmg**

It's best to create a new conda environment to use when running GMG. This will avoid any potential conflicts. e.g.:

    conda create -n gmgpy

Where the -n flag dictates what you want to name the new environment (you can call the environment whatever you like
but the documentation will use the name gmg-env).

On macOS and Linux this environment can then be activated using:

    conda activate gmgpy

Or on Windows use::

    activate gmgpy

Once the environment is activated any call to python will only "see" the packages install within the gmg-env
environment. To deactivate the environment either close the terminal window or on macOS and Linux run:

    conda deactivate

Or on Windows use:

    deactivate

**OS Support:**

GMG is currently only tested on macOS Monterey, Ubuntu 22.04 Linux and Windows 11 operating systems. It is possible
that issues may arise when trying to install on other distributions. Please raise any information regarding
installation problems on github.

Installing gmg
--------------

**NB: conda-forge and pip packages are currently in development and not yet deployed**

The gmg package is hosted on both conda-forge and PyPi (pip). It is highly recommended that
you install gmg in a standalone conda environment using option 1 below. This is ensure
you do not create any dependencies issues that may arise if you install system-wide using pip.

First, ensure your gmgpy conda environment is active. On macOS or Linux use:

    source activate gmgpy

Or on Windows use:

    activate gmgpy

**Installation option 1: using conda (highly recommended)**

    conda install gmgpy

**Installation option 2: using pip**

    pip install gmgpy

Launching gmg
-------------

after installation, make sure you are using a terminal with your gmgpy environment active 
(i.e., conda activate gmgpy), then, gmg can be launched from your terminal command line using:

    gmgpy
 
The splash screen and gmg model shell should appear. You can now begin using the software.

Tutorial and demo data
----------------------

Gmg ships with a tutorial dataset and benchmark demo models. If you installed gmg using the recommended miniforge 
approach as outlined above, then the gmg package directory is likely located in a path similar to:
    <pre>
    /Users/<strong>username</strong>/miniforge3/envs/gmgpy/lib/python3.<strong>X</strong>/site-packages/gmgpy
    </pre>

Within this directory are the two subdirectories demos/ and tutorial/. These can be copied or moved to another
location on your system or accessed from gmgpy directory when required. 

Getting started
---------------

The best way to get started is to read through the [Documentation](https://btozer.github.io/gmg/) starting with the [Manual](https://btozer.github.io/gmg/html/manual.html) and then work your way through the [Tutorial](https://btozer.github.io/gmg/html/tutorial.html). 
These can also be found by launching gmg and then navigating to:

    Help -> Documentation

Contributing
------------

Any contribution, big or small is welcome. Examples include:

1. Writing new code/functions.
2. Improving the user interface.
3. Submitting bug reports and spelling corrections.

They are all helpful!

**Becoming a developer**

The best way to become a developer is to create a folk of the github repository. You can then send a pull request once you are satisfied you want to share your modifications.

**Suggestions & bug reporting**

The best way to make a suggestion or open a bug report is to raise an issue on Github.

Please use the subject line:
    
    GMG suggestion: "short suggestion description"
    GMG bug report: "short bug description"

Please include as much detail as possible when describing the bug.

Citing
------

If you use GMG, please cite this github repository.
