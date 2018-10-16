---
title: 'GMG: An open-source two-dimensional geophysical modelling GUI'
tags:
  - Python
  - geophysics
  - forward modelling
  - gravity
  - magnetics
  - reflection seismic
authors:
  - name: Brook Tozer
    orcid: 0000-0003-0043-0580
    affiliation: 1
affiliations:
  - name: Scripps Institution of Oceanography, La Jolla, California 92093, USA
  - index: 1
date: 22 May 2018
bibliography: gmg.bib
---

# Summary
 
For decades, forward modelling of potential field data, such as gravity and magnetic 
anomalies, has been common practice within the geophysics community as a means for 
constraining subsurface structure. Many software packages (both freely available 
and commercially licenced) exist for performing such modelling. However, most, if not 
all of these packages suffer from at least one major drawback. Such draw-backs
include: (1) Being closed-source, with limited information regarding how anomalies are 
calculated; (2) the inability to calculate of both gravity and magnetic anomalies 
simultaneously; (3) providing no means for integrating complementary data within the 
modelling environment (e.g., earthquake hypocenters) and (4) being programmed in such a 
way that the software is cumbersome for integrating within a modern academic research 
project due to, for example, being a single platform release (usually Windows only), 
having poor I/O functionality and/or limited documentation.

``GMG`` is a Python package primarily intended as an interactive two-dimensional (2D) forward 
modelling GUI designed to resolve all of the drawbacks listed above. Both gravity and magnetic 
anomalies can be computed along a 2D profile consisting of subsurface bodies defined as any number 
of 2D polygons. Moreover, functions for displaying complementary data within the modelling environment, 
such as exploration well logs, seismic data and XY point data are provided. ``GMG`` has been designed 
with a minimalist user-interface and simple I/O functionality in order to enhance usability. The software 
is expected to be useful for both research purposes and for teaching of applied geophysics. 
Most importantly, ``GMG`` is open source, providing an environment where
users can add new functionality and optimise processes. In this way, it is hoped the 
software will naturally become more useful and streamlined over time.  

``GMG`` makes extensive use of functions from
the Scientic Computing in Python (SciPy) package [@oliphant2007]. In
particular, NumPy [@vanderwalt2011] data structures are used for data
handling and computational efficiency. Matplotlib [@hunter2007] plotting tools
are employed for displaying and interacting with graphics. The GUI is
implemented using the wxWidgets GUI toolkit, wxPython [@rappin2006].
Further software dependencies include Fatiando a Terra [@Uieda2013],
from which, the function fatiando.polygon is used for handling model layers and
ObsPy [@beyreuther2010], from which, the seismic plotting function
obspy.read is used for loading seismic data. The algorithms for calculating the
gravity and magnetic anomalies are adapted from [@bott1969] and [@talwani1964] respectively. 
The source code for ``GMG`` is stored on github at: https://github.com/btozer/gmg

# Acknowledgements

GMG was conceived at the University of Oxford and we thank Brook Keats for helpful discussion and 
contributions during development.

# References
