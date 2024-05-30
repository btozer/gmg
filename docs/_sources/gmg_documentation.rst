.. image:: _static/gmg_logo.png

.. raw:: html

    <div class="banner">
        <a href="http://www.github.com/btozer/gmg">Visit the GMG github repository</a>
    </div>

Overview
--------

GMG is an open-source Graphical User Interface (GUI) written in python 3, leveraging wxPython for GUI implementation. 

GMG is designed primarily for the integrated geophysical modelling of 2D potential field (gravity, vertical gravity 
gradient and magnetic) anomaly data. The package also allows the user to load and display seismic reflection data, 
well horizon tops, geological data and XY data (e.g., earthquake hypocenters).

GMG makes use of several other open-source python packages to preform various tasks including:

* WxPython
* Matplotlib
* SciPy
* NumPy
* ObsPy

(see `References <references.html>`_ for a detailed listing).

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

* Data I/O is designed to be as simple as possible using space delimited ASCII text files.
* Import and display observed topography, gravity, vertical gravity gradient and magnetic data.
* Ability to apply filters to observed data within the modelling environment.
* Import and display seismic reflection data.
* Add and manipulate model layers (subsurface bodies) using a simple, responsive interactive interface.
* “Pinch out/snap” layers against adjacent layers.
* Calculate the predicted potential field anomalies produced by any combination of model layers (i.e., quickly include/exclude selected layers). 
* Display well horizon tops.
* Display XY data (e.g., earthquake hypocenters or geological surface contacts).
* Export model data (e.g., predicted anomalies and layer geometries) as ASCII text files.
* Save model figures as vector graphics in Portable Document Format (PDF) file.
* Access to an on-the-fly python terminal for access to all underlying python objects.
* **The ability to modify the software in any way the user wishes.**

**NB: GMG is currently in development and a work in progress. As such, some documentation may be 
incomplete and some features may not work as expected. If you have issues, please reach 
out so these can be corrected or the documentation modified.**

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Getting Started

    installation
    getting_started

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: User Guide

    manual
    tutorial
    contribute

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Credits

    references
    related_software
    licence

Release Notes
-------------

Version: |version|

Last updated: |today|

Notices:
--------

If you would like to make a comment, create a feature request or report a bug please use the github page.

Citing:
-------

If you use GMG, please cite this github repository (Citation to come once the package is further developed).
