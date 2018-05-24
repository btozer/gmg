#!/bin/bash

pandoc paper.md --latex-engine=xelatex --bibliography paper.bib -o paper.pdf
