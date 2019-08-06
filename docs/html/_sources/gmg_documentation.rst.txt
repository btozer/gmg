.. image:: _static/gmg_logo.png

.. raw:: html

    <div class="banner">
        <a href="http://www.github.com/btozer/gmg">Visit the GMG github repository</a>
    </div>

Overview
--------

GMG is an open-source Graphical User Interface (GUI) designed principally for modelling
2D potential field (gravity and magnetic) profiles. GMG is written in python 3.7 and uses
wxPython for GUI implementation.

GMG makes use of several other open-source python packages to preform various tasks
(see `References <references.html>`_ for a detailed listing).

The project was instigated after failing to find an adequate open-source package (in which
the source code can be viewed and modified by the user) for performing 2D geophysical
modeling tasks. Inspiration came from fatiando a terra and GMT.

It is anticipated that GMG will also be useful for teaching purposes.

**NB: GMG is currently in development. Some documentation is incomplete and some features may not work as expected.**

Key features
------------
* Data I/O is made as simple as possible using space delimited ASCII text files.
* Import and display observed topography, gravity and magnetic data.
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
* Save model figures as vector graphics in Portable Document Format (PDF) file.

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
Alternatively, you can send an e-mail to: btozer@ucsd.edu

Citing:
-------

If you use GMG, please cite this github repository (Citation to come once the package is further developed).
