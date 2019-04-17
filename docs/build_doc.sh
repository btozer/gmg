#!/bin/bash

#% HTML
rm -r  html/
sphinx-build -b html _sources/ html/

#% LATEX
#rm -r  _build/latex
#sphinx-build -b latex _sources/ _build/latex
