#!/bin/bash

## IF INPUT ARG IS NOT PROVIDED
if [ -z "$1" ]; then
    echo "Please provide action as input an argument"
    echo "Usage: ./build_pip.sh [build|submit_test|submit_final|clean]"
    exit 1
fi

## BUILD
if [ "$1" == "build" ]; then
    rm -r dist/
    python3 -m build
    pip install .
    exit 0
fi

## SUBMIT TEST
if [ "$1" == "submit_test" ]; then
    python3 -m twine upload --repository testpypi dist/*
    exit 0
fi

## SUBMIT FINAL
if [ "$1" == "submit_final" ]; then
    python3 -m twine upload --repository pypi dist/*
    exit 0
fi

## CLEAN
if [ "$1" == "clean" ]; then
    rm -r dist/
    rm -r build/
    rm -r *.egg-info/
    exit 0
fi