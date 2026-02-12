import os
import numpy as np
import netCDF4
import xarray as xr
import pygrib
from osgeo import gdal
import glob
import random
import string
import sys
import gc
import psutil

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
def runner(maxlvl,insurface,outfile):
    try:
        if 'C1' in insurface:
            inprofile=insurface.replace('C1','C2')
        elif 'C3' in insurface:
            inprofile=insurface.replace('C3','C4')
            inprofile=inprofile[0:-3]+'*'
            inprofile=(glob.glob(inprofile))[0]
        else:
            print('Unrecognised surface parameter file: '+\
                  os.path.basename(insurface))
            exit()
    except:
        print('Failed to find profile file matching surface file?: '+insurface)
        return

    print('Opening file: '+insurface)
    grbs		=	pygrib.open(insurface)
    print('Opening file: '+inprofile)
    grbp        =   pygrib.open(inprofile)

    print('	-	Determining array sizes and allocating space.')
    if maxlvl < 0:
        for grb in grbs:
            if (grb.name!="Temperature"):	continue
            maxlvl+=1
        grbs.close()
        grbs	=	pygrib.open(insurface)

    # The following doesn't seem to work with the forecast 
    #files supplied be CEDA
    #lats		=	grbs[1].distinctLatitudes
    #lons		=	grbs[1].distinctLongitudes
    # Could use the built in informatioin from the file to
    # build the grids, but we don't need to - all we need are
    # the dimensions
    # lats,lons = grbs[0].latlons()

    #Varnames:
    #	Q,	T,	O3
    #	LNSP,Z,	CI,	ASN,	TCWV,	SD,	U10,	V10,	T2,	SKT,	SSTK,	AL
    shaper	=	grbs[1].values.shape

    # Construct our own lat/lon arrays, based on the number of elements
    # in the data and what we know of the way ECMWF store the coordinates
    # lats = (np.arange(shaper[0]/2, -1*shaper[0]/2, -1) - 0.5)*180.0/shaper[0] 
    # lons = np.arange(0, shaper[1], 1) * 360.0/shaper[1]
    lats,lons = grbs[1].latlons()
    lats = lats[:,0]
    lons = lons[0,:]
    #print 'Latitude array has dimension: ',lats.__len__()
    #print lats[0], ',',lats[1], ',', lats[2], '...',\
    #      lats[lats.__len__()-3], ',', lats[lats.__len__()-2], ',',\
    #      lats[lats.__len__()-1]
    #print 'Longitude array has dimension: ',lons.__len__()
    #print lons[0], ',',lons[1], ',', lons[2], '...',\
    #      lons[lons.__len__()-3], ',', lons[lons.__len__()-2], ',',\
    #      lons[lons.__len__()-1]

    # Profile variables
    q		=	np.zeros((maxlvl,shaper[0],shaper[1]))
    t		=	np.zeros((maxlvl,shaper[0],shaper[1]))
    o3		=	np.zeros((maxlvl,shaper[0],shaper[1]))
    # Surface variables
    lnsp	=	np.zeros((shaper[0],shaper[1]))
    z		=	np.zeros((shaper[0],shaper[1]))
    ci		=	np.zeros((shaper[0],shaper[1]))
    asn		=	np.zeros((shaper[0],shaper[1]))
    tcwv	=	np.zeros((shaper[0],shaper[1]))
    sd		=	np.zeros((shaper[0],shaper[1]))
    u10		=	np.zeros((shaper[0],shaper[1]))
    v10		=	np.zeros((shaper[0],shaper[1]))
    t2		=	np.zeros((shaper[0],shaper[1]))
    skt		=	np.zeros((shaper[0],shaper[1]))
    sstk	=	np.zeros((shaper[0],shaper[1]))
    al		=	np.zeros((shaper[0],shaper[1]))


    print('	-	Retrieving data from the profile GRIB file.')
    for grb in grbp:
        if (grb.name=="Temperature"):
            t[grb.level-1,:,:]=grb.values
            continue
        if (grb.name=="Specific humidity"):
            q[grb.level-1,:,:]=grb.values
            continue
        if (grb.name=="Ozone mass mixing ratio"):
            o3[grb.level-1,:,:]=grb.values
            continue
    print('	-	Retrieving data from the surface GRIB file.')
    for grb in grbs:
        if (grb.name=="Logarithm of surface pressure"):
            lnsp[:,:]=grb.values
            pts=(lnsp==grb.missingValue).nonzero()
            lnsp[pts]=grb['minimum']
            continue
        if (grb.name=="Surface pressure"):
            lnsp[:,:]=np.log(grb.values)
            pts=(lnsp==np.log(grb.missingValue)).nonzero()
            lnsp[pts]=grb['minimum']
            continue
        if (grb.name=="Geopotential"):
            z[:,:]=grb.values
            continue
        if (grb.name=="Sea ice area fraction"): # was Sea-ice cover
            ci[:,:]=grb.values
            pts=(ci==grb.missingValue).nonzero()
            ci[pts]=grb['minimum']
            continue
        if (grb.name=="Snow albedo"):
            asn[:,:]=grb.values
            pts=(asn==grb.missingValue).nonzero()
            asn[pts]=grb['minimum']
            continue
        if (grb.name=="Total column water vapour"):
            tcwv[:,:]=grb.values
            pts=(tcwv==grb.missingValue).nonzero()
            tcwv[pts]=grb['minimum']
            continue
        if (grb.name=="Snow depth"):
            sd[:,:]=grb.values
            pts=(sd==grb.missingValue).nonzero()
            sd[pts]=grb['minimum']
            continue
        if (grb.name=="10 metre U wind component"):
            u10[:,:]=grb.values
            pts=(u10==grb.missingValue).nonzero()
            u10[pts]=grb['minimum']
            continue
        if (grb.name=="10 metre V wind component"):
            v10[:,:]=grb.values
            pts=(v10==grb.missingValue).nonzero()
            v10[pts]=grb['minimum']
            continue
        if (grb.name=="2 metre temperature"):
            t2[:,:]=grb.values
            pts=(t2==grb.missingValue).nonzero()
            t2[pts]=grb['minimum']
            continue
        if (grb.name=="Skin temperature"):
            skt[:,:]=grb.values
            pts=(skt==grb.missingValue).nonzero()
            skt[pts]=grb['minimum']
            continue
        if (grb.name=="Sea surface temperature"):
            sstk[:,:]=grb.values
            pts=(sstk==grb.missingValue).nonzero()
            sstk[pts]=grb['minimum']
            continue
        if (grb.name=="Land-sea mask"):
            al[:,:]=grb.values
            pts=(al==grb.missingValue).nonzero()
            al[pts]=grb['minimum']
            continue

    print('	-	Saving variables into ORAC-style NetCDF file.')
    ncfile	=	netCDF4.Dataset(outfile, 'w',format='NETCDF3_64BIT')
    ncfile.createDimension('latitude',shaper[0])
    ncfile.createDimension('longitude',shaper[1])
    ncfile.createDimension('hybrid',maxlvl)
    lato	=	ncfile.createVariable('latitude','f',('latitude'),
                                      fill_value=-999)
    lato[:]	=	lats
    lono	=	ncfile.createVariable('longitude','f',('longitude'),
                                      fill_value=-999)
    lono[:]	=	lons
    hybo	=	ncfile.createVariable('hybrid','f',('hybrid'),
                                      fill_value=-999)
    hybo[:]	=	np.arange(maxlvl)

    qv		=	ncfile.createVariable('Q','f',('hybrid','latitude',
                                      'longitude'),fill_value=2.0e20)
    qv[:]	=	q
    tv		=	ncfile.createVariable('T','f',('hybrid','latitude',
                                      'longitude'),fill_value=2.0e20)
    tv[:]	=	t
    o3v		=	ncfile.createVariable('O3','f',('hybrid','latitude',
                                          'longitude'),fill_value=2.0e20)
    o3v[:]	=	o3

    zv		=	ncfile.createVariable('Z','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    zv[:]	=	z
    lnspv	=	ncfile.createVariable('LNSP','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    lnspv[:]=	lnsp
    civ		=	ncfile.createVariable('CI','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    civ[:]	=	ci
    asnv	=	ncfile.createVariable('ASN','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    asnv[:]	=	asn
    tcwvv	=	ncfile.createVariable('TCWV','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    tcwvv[:]=	tcwv
    sdv		=	ncfile.createVariable('SD','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    sdv[:]	=	sd
    u10v	=	ncfile.createVariable('U10','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    u10v[:]	=	u10
    v10v	=	ncfile.createVariable('V10','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    v10v[:]	=	v10
    t2v		=	ncfile.createVariable('T2','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    t2v[:]	=	t2
    sktv	=	ncfile.createVariable('SKT','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    sktv[:]	=	skt
    sstkv	=	ncfile.createVariable('SSTK','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    sstkv[:]=	sstk
    alv		=	ncfile.createVariable('AL','f',('latitude','longitude'),
                                      fill_value=2.0e20)
    alv[:]	=	al
    
    grbs.close()

    return

def proc_dir(indir, outdir, infile=None):
    infiles = glob.glob(indir + '/C*')
    infiles.sort()

    files = []
    if infile:
        f2 = os.path.splitext(os.path.basename(infile))[0]
        files.append(f2)
    else:
        for f in infiles:
            f2 = os.path.splitext(os.path.basename(f))[0]
            files.append(f2)

    for f in files:
        proc = psutil.Process(os.getpid())
        print(f"Open files: {len(proc.open_files())}", proc.open_files())
        final_output_file = os.path.join(outdir, f + '.nc')
        
        if not os.path.isfile(final_output_file):
            # Generate a temporary output file with a random 4-letter suffix
            random_suffix = ''.join(random.choices(string.ascii_lowercase, k=4))
            temp_output_file = final_output_file + f".{random_suffix}"
            
            try:
                # Call runner with the temporary output file
                runner(137, os.path.join(indir, f), temp_output_file)
                print(temp_output_file)
                # Rename temporary file to final file after processing
                os.rename(temp_output_file, final_output_file)
                print(f"Renamed: {temp_output_file} -> {final_output_file}")
            except Exception as e:
                print(f"Error processing {f}: {e}")
                # Clean up temporary file in case of failure
                if os.path.exists(temp_output_file):
                    os.remove(temp_output_file)
            finally:
                gc.collect()
    return 0

if __name__=='__main__':
    indir  = sys.argv[1]
    outdir = sys.argv[2]
    if len(sys.argv) > 3:
        s = proc_dir(indir,outdir,sys.argv[3])
    else:
        s = proc_dir(indir,outdir)
    sys.exit(s)
