#!/bin/bash

pandoc paper.md --pdf-engine=xelatex --bibliography paper.bib -o paper.pdf
