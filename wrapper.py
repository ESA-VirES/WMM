#!/usr/bin/env python
import sys
import argparse
import os
import os.path
import subprocess
import numpy
import gdal
import osr
import ogr

def main(args):

    parser = argparse.ArgumentParser()
    parser.add_argument("resolution", nargs=1, type=float)
    parser.add_argument("height", nargs=1, type=int)
    parser.add_argument("time", nargs=1, type=float)
    parser.add_argument("product", nargs=1, type=int)
    parser.add_argument("interval", nargs=1, type=int)

    parsed = parser.parse_args(args)

    exe='./wmm_grid.exe'
    proc = subprocess.Popen([exe],1,exe,subprocess.PIPE,sys.stdout,subprocess.STDOUT) 

    # Minimum Latitude (in decimal degrees)
    latmin = -90

    # Maximum Latitude (in decimal degrees)
    latmax = 90

    # Minimum Longitude (in decimal degrees)
    lonmin = -180

    # Maximum Longitude (in decimal degrees)
    lonmax = 180

    # Step Size (in decimal degrees)
    deg_interval = parsed.resolution[0]

    # 1: above mean sea level; 2: WGS-84 ellipsoid
    height = "2"

    # Minimum Height above the WGS-84 Ellipsoid (in km)
    heightmax = parsed.height[0]

    # Maximum Height above the WGS-84 Ellipsoid (in km)
    heightmin = parsed.height[0]

    # height step size (in km)
    height_interval = 0

    # decimal year starting time
    timestart = parsed.time[0]

    # decimal year ending time
    timeend = parsed.time[0]

    # time step size
    time_interval = "0"

    # geomagnetic element to print. Your options are :
    # 1. Declination	9.   Ddot
    # 2. Inclination	10. Idot
    # 3. F				11. Fdot
    # 4. H				12. Hdot
    # 5. X				13. Xdot
    # 6. Y				14. Ydot
    # 7. Z				15. Zdot    
    # 8. GV				16. GVdot
    product = parsed.product[0]

    # select output (1 for file)
    output = "1"

    # output filename
    tempfile = "temp_result"

    # generate geotransform values
    resolution = deg_interval
    rasterxsize = (lonmax-lonmin)/deg_interval
    rasterysize = (latmax-latmin)/deg_interval

    geotransform = [
        lonmin,
        resolution,
        0,
        latmax,
        0,
        -resolution
    ]

    wmmxmin = latmin+deg_interval/2
    wmmxmax = latmax-deg_interval/2
    wmmymin = lonmin+deg_interval/2
    wmmymax = lonmax-deg_interval/2

    print >>proc.stdin, wmmxmin
    print >>proc.stdin, wmmxmax
    print >>proc.stdin, wmmymin
    print >>proc.stdin, wmmymax
    print >>proc.stdin, deg_interval
    print >>proc.stdin, height
    print >>proc.stdin, heightmax
    print >>proc.stdin, heightmin
    print >>proc.stdin, height_interval
    print >>proc.stdin, timestart
    print >>proc.stdin, timeend
    print >>proc.stdin, time_interval
    print >>proc.stdin, product
    print >>proc.stdin, output
    print >>proc.stdin, tempfile
    print >>proc.stdin

    status = proc.wait()
    print "STATUS:", status 
    proc.stdin.close()

    #print wmmxmin, wmmxmax, wmmymin, wmmymax

    # create 1d numpy array
    raster_out = numpy.zeros((rasterxsize*rasterysize), "float32")
    linenumber = 0

    with open(tempfile) as f:
        for line in f:
            line_array = line.split( )
            raster_out[linenumber] = str(line_array[4])
            linenumber = linenumber+1

    raster_2d = numpy.flipud(raster_out.reshape(rasterysize,rasterxsize))

    #print raster_2d.shape
    #print rasterxsize, rasterysize

    driver = gdal.GetDriverByName('MEM')
    mem_ds = driver.Create('', int(rasterxsize), int(rasterysize), 1, gdal.GDT_Float32)
    mem_ds.GetRasterBand(1).WriteArray(raster_2d)
    mem_ds.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    mem_ds.SetProjection( srs.ExportToWkt() )

    contour_steps = parsed.interval[0]

    # clean up from previous runs
    try:
        os.remove('contour.shp')
    except:
        pass
    try:
        os.remove('contour.dbf')
    except:
        pass
    try:
        os.remove('contour.shx')
    except:
        pass

    ogr_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource('contour.shp')
    ogr_lyr = ogr_ds.CreateLayer('contour')
    field_defn = ogr.FieldDefn('ID', ogr.OFTInteger)
    ogr_lyr.CreateField(field_defn)
    field_defn = ogr.FieldDefn('elev', ogr.OFTReal)
    ogr_lyr.CreateField(field_defn)

    gdal.ContourGenerate(mem_ds.GetRasterBand(1), contour_steps, 0, [], 0, 0, ogr_lyr, 0, 1)

    os.remove(tempfile)

if __name__ == "__main__":
        main(sys.argv[1:])