#!/bin/bash
in=$1
out=$2
size=$3

convert ${in}.png -resize ${size}x${size} ${out}.png
