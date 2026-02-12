#!/usr/bin/env python3
import xarray as xr
import datetime
from getopt import getopt
import sys
import os

final_vars = [
'time',
'lat',
'lon',
'boa_swdn_tot',
'boa_par_tot',
'boa_par_dif',
'boa_swdn_dif',
'stemp']

options, operands = getopt(sys.argv[1:], "", ["msi_root=", 'bugsradfile='])
for o,v in options:
    if o == "--msi_root":
        msi_root = v
    elif o == "--bugsradfile":
        bugsradfile = v

bugs = xr.open_dataset(bugsradfile, engine='netcdf4', decode_cf=False)
bugs.attrs['Project'] = 'SISEM'
bugs.attrs['File_Name'] = bugsradfile.replace('bugsrad', 'sisem')
bugs.attrs['Creator_url'] = 'https://github.com/ORAC-CC/orac'
bugs.attrs['License'] = 'https://github.com/ORAC-CC/orac/wiki/License'
bugs.attrs['title'] = 'SISEM L2 Bugsrad File'
bugs[final_vars].to_netcdf(bugsradfile.replace('bugsrad', 'sisem'))
