
#GMG: An open source two-dimensional geophysical modelling GUI

---
tags:
  - Python
  - geophysics
  - forward modelling
  - gravity
  - magnetics
  - reflection seismic
author:
  - name: Brook Tozer
    orcid: 0000-0003-0043-0580
affiliation:
 - name: Scripps Institution of Oceanography, La Jolla, California 92093, USA
---

# Summary
For decades, forward modelling of potential field data, such as 
gravity and magnetic anomalies, has been common practice within the 
geophysics community as a means of constraining subsurface 
structure. Many software packages (both freely available and commercially licenced) 
exist for preforming such modelling. However, most, if not
all of these packages, suffer from at least one major draw-back. Such draw-backs
include: (1) being closed source; (2) not allowing for the calculation of both
gravity and magnetic anomalies; (3) Not facilitating the integration of
complimentary data within the modelling environment and (4) being programmed in 
such a way that the software is  cumbersome for integrating within an 
academic research project due to, for example, (5) being a single platform 
release (usually Windows only), (6) having poor I/O functionality and (7) poor 
documentation.

``GMG`` is an open-source Python package primarily intended as an
interactive, user-friendly two-dimensional geophysical forward modelling GUI
that resolves all of the draw-backs listed above. Both gravity and magnetic
anomalies can be computed along a 2D profile consisting of subsurface bodies
defined as any number of 2D polygons. Moreover, functions for displaying
complementary data within the modelling environment, such as exploration well
logs and seismic data, are provided. ``GMG`` has been designed with a minimalist
user-interface and simple I/O in order to enhance usability and it is expected
to be useful for both researchers and teaching of exploration geophysics. Most
importantly ``GMG`` is fully open source. Hence, users can add new 
functionality that is beneficial to the wider community and the software should
naturally become more useful and streamlined over time.  

``GMG`` makes extensive use of functions from
the Scientic Computing in Python (SciPy) package (Oliphant, 2007). In
particular, NumPy (Van Der Walt et al., 2011) data structures are used for data
handling and computational efficiency. Matplotlib (Hunter, 2007) plotting tools
are employed for displaying and interacting with graphics. The GUI is
implemented using the wxWidgets GUI toolkit, wxPython (Rappin and Dunn, 2006).
Further, software dependencies include Fatiando a Terra (Uieda et al., 2013),
from which, the function fatiando.polygon is used for handling model layers and
ObsPy (Beyreuther et al., 2010), from which, the seismic plotting function
obspy.read is used for loading seismic data. The source code for ``GMG`` is
stored on github at: https://github.com/btozer/gmg

# Acknowledgements

GMG was conceived at the University of Oxford and I acknowledge contributions
from Brook Keats.

# References
bibliography: bibliography.bib



Figures can be included like this: ![Example figure.](figure.png)

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.