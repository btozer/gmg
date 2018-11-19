#!/bin/bash

#% HTML
rm -r  _build/html
sphinx-build -b html _sources/ ./

#% LATEX
#rm -r  _build/latex
#sphinx-build -b latex _sources/ _build/latex
