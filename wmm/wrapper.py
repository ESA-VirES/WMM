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
from uuid import uuid4
import logging


logger = logging.getLogger(__name__)

zoomresolution = {
    0.70312500:0,
    0.35156250:1,
    0.17578125:2,
    0.08789063:3,
    0.04394531:4,
    0.02197266:5,
    0.01098633:6,
    0.00549316:7,
    0.00274658:8,
    0.00137329:9,
    0.00068665:10
}

# product names:(id, iso_interval flavor)
products = {
    "Decl":(1, 1),
    "Incl":(2, 1),
    "F":(3, 2),
    "H":(4, 2),
    "X":(5, 2),
    "Y":(6, 2),
    "Z":(7, 2)
    #"GV":{8,
    #"Ddot":{9.,
    #"Idot":{10,
    #"Fdot":{11,
    #"Hdot":{12,
    #"Xdot":{13,
    #"Ydot":{14,
    #"Zdot":{15,
    #"GVdot":{16
}

iso_intervals = {
    1:{
        0: 10,
        1: 10,
        2: 5,
        3: 5,
        4: 1,
        5: 1,
        6: 0.5,
        7: 0.5,
        8: 0.1
    },
    2:{
        0: 5000,
        1: 5000,
        2: 2500,
        3: 1000,
        4: 1000,
        5: 1000,
        6: 500,
        7: 500,
        8: 100
    }
}


def main(filename, product, bbox, pixelsize, height, time):
    # configuration
    ## resolutions (in degrees) per zoom level for EPSG 4326
    
    exe = os.path.join(os.path.dirname(__file__), 'wmm_grid.exe')
    
    """if not os.path.exists("WMM.COF"):
        raise Exception(
            "Coefficient file missing: %s" 
            % os.path.join(os.path.dirname(__file__), "WMM.COF")
        )"""

    proc = subprocess.Popen(
        [exe], 1, exe, subprocess.PIPE, open(os.devnull),
        subprocess.STDOUT, cwd=os.path.dirname(__file__)
    )

    latmin, lonmin, latmax, lonmax = bbox

    # Step Size (in decimal degrees)
    #deg_interval = resolution
    deg_interval = (latmax-latmin)/pixelsize

    # Minimum Height above the WGS-84 Ellipsoid (in km)
    heightmax = height

    # Maximum Height above the WGS-84 Ellipsoid (in km)
    heightmin = height

    # height step size (in km)
    height_interval = 0

    # decimal year starting time
    timestart = time

    # decimal year ending time
    timeend = time

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
    product_id = products[product][0]

    # select output (1 for file)
    output = "1"

    # output filename
    tempfile = "/tmp/%s" % uuid4().hex[:10]

    # generate geotransform values
    geotransform = [
        lonmin,
        deg_interval,
        0,
        latmax,
        0,
        -deg_interval
    ]

    wmmxmin = lonmin-deg_interval/2
    wmmxmax = lonmax+deg_interval/2
    wmmymin = latmin-deg_interval/2
    wmmymax = latmax+deg_interval/2

    print >>proc.stdin, wmmymin
    print >>proc.stdin, wmmymax
    print >>proc.stdin, wmmxmin
    print >>proc.stdin, wmmxmax
    print >>proc.stdin, deg_interval
    print >>proc.stdin, height
    print >>proc.stdin, heightmax
    print >>proc.stdin, heightmin
    print >>proc.stdin, height_interval
    print >>proc.stdin, timestart
    print >>proc.stdin, timeend
    print >>proc.stdin, time_interval
    print >>proc.stdin, product_id
    print >>proc.stdin, output
    print >>proc.stdin, tempfile
    print >>proc.stdin

    status = proc.wait()
    print "STATUS:", status 
    proc.stdin.close()

    values = []
    xvalues = set()
    yvalues = set()

    with open(tempfile) as f:
        for line in f:
            line_array = line.split( )
            values.append(str(line_array[4]))
            xvalues.add(str(line_array[1]))
            yvalues.add(str(line_array[0]))

    rasterxsize = len(xvalues)
    rasterysize = len(yvalues)

    logger.info("size %d %d " %(rasterxsize, rasterysize))

    os.remove(tempfile)

    # create 1d numpy array, reshape to 2d, and flip to correct order
    raster_out = numpy.fromiter(values, "float32")
    raster_2d = numpy.flipud(raster_out.reshape(rasterysize,rasterxsize))

    #print raster_2d.shape
    #print rasterxsize, rasterysize

    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(
        filename, rasterxsize, rasterysize, 1, gdal.GDT_Float32
    )
    ds.GetRasterBand(1).WriteArray(raster_2d)
    ds.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())

    return ds

    #contour_steps = parsed.interval[0]
    """
    zoom = round(deg_interval, 8)
    contour_steps = iso_intervals[
        products[product][1]
    ][zoomresolution[zoom]]

    #print "zoomlevel: " + zoomresolution[deg_interval]
    #print iso_intervals[products[parsed.product[0]][1]][zoomresolution[deg_interval]]

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


    ogr_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(filename)
    ogr_lyr = ogr_ds.CreateLayer('contour')
    field_defn = ogr.FieldDefn('ID', ogr.OFTInteger)
    ogr_lyr.CreateField(field_defn)
    field_defn = ogr.FieldDefn('elev', ogr.OFTReal)
    ogr_lyr.CreateField(field_defn)

    gdal.ContourGenerate(mem_ds.GetRasterBand(1), contour_steps, 0, [], 0, 0, ogr_lyr, 0, 1)

    #print min(raster_out), max(raster_out)

    """


if __name__ == "__main__":
    args = sys.argv[1:]
    parser = argparse.ArgumentParser()
    parser.add_argument("bbox", nargs=4, type=float)
    #parser.add_argument("resolution", nargs=1, type=float)
    parser.add_argument("pixelsize", nargs=1, type=int)
    parser.add_argument("height", nargs=1, type=int)
    parser.add_argument("time", nargs=1, type=float)
    parser.add_argument("product", nargs=1, type=str)
    parser.add_argument("--contour", action="store_true", default=False)
    #parser.add_argument("interval", nargs=1, type=float)

    parsed = parser.parse_args(args)

    main(
        "contour.shp", parsed.product[0], parsed.bbox, parsed.pixelsize[0], 
        parsed.height[0], parsed.time[0], parsed.contour
    )
