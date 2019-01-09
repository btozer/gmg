"""
Cut 2D profile from SEAM Interpretation Challenge 3D volume

https://seg.org/News-Resources/Research-and-Data/SEAM

https://wiki.seg.org/wiki/SEAM
"""

import sys
from obspy import read, Stream

## READ INPUT FILE NAME
input = sys.argv[1]

## READ IN SEGY DATA
section = read(input, unpack_trace_headers=True)

## EXTRACT 2D PROFILE
profile = Stream()
for i in range(0, len(section.traces)):
    if section.traces[i].stats['segy']['trace_header']['group_coordinate_y'] == 24000:
        print(i)
        profile += section.traces[i]

## WRITE 2D PROFILE TO DISC
profile.write("2d_profile_24000N.segy", format="SEGY")
