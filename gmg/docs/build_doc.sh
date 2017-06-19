#!/bin/bash

#% HTML
sphinx-build -b html _sources/ _build/html/

#% LATEX
sphinx-build -b latex _sources/ _build/latex
