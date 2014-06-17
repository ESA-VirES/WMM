from datetime import datetime, date, timedelta
import time
import numpy as np
cimport numpy as np

cdef extern from "../GeomagnetismHeader.h":
    ctypedef struct MAGtype_MagneticModel:
        double EditionDate
        double epoch
        char ModelName[32]
        double *Main_Field_Coeff_G
        double *Main_Field_Coeff_H
        double *Secular_Var_Coeff_G
        double *Secular_Var_Coeff_H
        int nMax
        int nMaxSecVar
        int SecularVariationUsed
        double CoefficientFileEndDate

    ctypedef struct MAGtype_Ellipsoid:
        pass

    ctypedef struct MAGtype_CoordGeodetic:
        double lambda_ "lambda"
        double phi
        double HeightAboveEllipsoid
        double HeightAboveGeoid
        int UseGeoid

    ctypedef struct MAGtype_CoordSpherical:
        double lambda_ "lambda"
        double phig
        double HeightAboveGeoid

    ctypedef struct MAGtype_Date:
        int Year
        int Month
        int Day
        double DecimalYear

    ctypedef struct MAGtype_LegendreFunction:
        pass

    ctypedef struct MAGtype_MagneticResults:
        pass

    ctypedef struct MAGtype_SphericalHarmonicVariables:
        pass

    ctypedef struct MAGtype_GeoMagneticElements:
        double Decl
        double Incl
        double F
        double H
        double X
        double Y
        double Z
        double GV
        double Decldot
        double Incldot
        double Fdot
        double Hdot
        double Xdot
        double Ydot
        double Zdot
        double GVdot

    ctypedef struct MAGtype_Geoid:
        int NumbGeoidCols
        int NumbGeoidRows
        int NumbHeaderItems
        int ScaleFactor
        float *GeoidHeightBuffer
        int NumbGeoidElevs
        int Geoid_Initialized
        int UseGeoid

    ctypedef struct MAGtype_CoordGeodeticStr:
        pass

    ctypedef struct MAGtype_UTMParameters:
        pass

    int MAG_SetDefaults(MAGtype_Ellipsoid *Ellip, MAGtype_Geoid *Geoid)
    int MAG_robustReadMagModels(char *filename, MAGtype_MagneticModel *(*magneticmodels)[], int array_size)
    MAGtype_MagneticModel *MAG_AllocateModelMemory(int NumTerms)
    MAGtype_SphericalHarmonicVariables *MAG_AllocateSphVarMemory(int nMax)
    MAGtype_LegendreFunction *MAG_AllocateLegendreFunctionMemory(int NumTerms)
    int MAG_ConvertGeoidToEllipsoidHeight(MAGtype_CoordGeodetic *CoordGeodetic, MAGtype_Geoid *Geoid)
    int MAG_GeodeticToSpherical(MAGtype_Ellipsoid Ellip, MAGtype_CoordGeodetic CoordGeodetic, MAGtype_CoordSpherical *CoordSpherical)
    int MAG_ComputeSphericalHarmonicVariables(MAGtype_Ellipsoid Ellip,
        MAGtype_CoordSpherical CoordSpherical,
        int nMax,
        MAGtype_SphericalHarmonicVariables * SphVariables)
    int MAG_AssociatedLegendreFunction(MAGtype_CoordSpherical CoordSpherical, int nMax, MAGtype_LegendreFunction *LegendreFunction)
    int MAG_TimelyModifyMagneticModel(MAGtype_Date UserDate, MAGtype_MagneticModel *MagneticModel, MAGtype_MagneticModel *TimedMagneticModel)
    int MAG_Summation(MAGtype_LegendreFunction *LegendreFunction,
        MAGtype_MagneticModel *MagneticModel,
        MAGtype_SphericalHarmonicVariables SphVariables,
        MAGtype_CoordSpherical CoordSpherical,
        MAGtype_MagneticResults *MagneticResults)
    int MAG_SecVarSummation(MAGtype_LegendreFunction *LegendreFunction,
        MAGtype_MagneticModel *MagneticModel,
        MAGtype_SphericalHarmonicVariables SphVariables,
        MAGtype_CoordSpherical CoordSpherical,
        MAGtype_MagneticResults *MagneticResults)
    int MAG_SecVarSummationSpecial(MAGtype_MagneticModel *MagneticModel,
        MAGtype_SphericalHarmonicVariables SphVariables,
        MAGtype_CoordSpherical CoordSpherical,
        MAGtype_MagneticResults *MagneticResults)
    int MAG_RotateMagneticVector(MAGtype_CoordSpherical,
        MAGtype_CoordGeodetic CoordGeodetic,
        MAGtype_MagneticResults MagneticResultsSph,
        MAGtype_MagneticResults *MagneticResultsGeo)
    int MAG_CalculateGeoMagneticElements(MAGtype_MagneticResults *MagneticResultsGeo, MAGtype_GeoMagneticElements *GeoMagneticElements)
    int MAG_CalculateSecularVariationElements(MAGtype_MagneticResults MagneticVariation, MAGtype_GeoMagneticElements *MagneticElements)
    int MAG_FreeLegendreMemory(MAGtype_LegendreFunction *LegendreFunction)
    int MAG_FreeMagneticModelMemory(MAGtype_MagneticModel *MagneticModel)
    int MAG_FreeSphVarMemory(MAGtype_SphericalHarmonicVariables *SphVar)



cdef extern from "../EGM9615.h":
    float GeoidHeightBuffer[]


"""
cdef class MagneticModel:
    cdef MAGtype_MagneticModel * _magnetic_model[1]
    cdef MAGtype_Ellipsoid _ellipsoid
    cdef MAGtype_Geoid _geoid

    def __cinit__(self, char* filename="WMM.COF"):
        #MAG_robustReadMagModels(filename, &self._magnetic_model, 1)

        MAG_SetDefaults(&self._ellipsoid, &self._geoid)
"""


cdef class Sampler:
    cdef MAGtype_MagneticModel * _magnetic_model[1]
    cdef MAGtype_Ellipsoid _ellipsoid
    cdef MAGtype_Geoid _geoid

    def __cinit__(self, char* filename="WMM.COF"):
        #cdef MAGtype_MagneticModel * magnetic_model[1]# =  self._magnetic_model
        
        MAG_robustReadMagModels(filename, <MAGtype_MagneticModel *(*)[]>&(self._magnetic_model[0]), 1)
        MAG_SetDefaults(&self._ellipsoid, &self._geoid)

    def __dealloc__(self):
        pass


PRODUCTS = [
    "Decl",
    "Incl",
    "F",
    "H",
    "X",
    "Y",
    "Z",
    "GV",
    "Decldot",
    "Incldot",
    "Fdot",
    "Hdot",
    "Xdot",
    "Ydot",
    "Zdot",
    "GVdot"
]

def grid(double minx, double miny, double minz, double mint,
         double maxx, double maxy, double maxz, double maxt,
         double dx,   double dy,   double dz,   double dt, 
         char* product, const char* filename
         ):
    cdef MAGtype_MagneticModel * magnetic_models[1]
    cdef MAGtype_MagneticModel * magnetic_model
    cdef MAGtype_Ellipsoid ellipsoid
    cdef MAGtype_Geoid geoid

    cdef int product_id = PRODUCTS.index(product) + 1
    cdef double value = 0.0


    MAG_robustReadMagModels(filename, &magnetic_models, 1)
    magnetic_model = magnetic_models[0]
    MAG_SetDefaults(&ellipsoid, &geoid)

    geoid.GeoidHeightBuffer = GeoidHeightBuffer;
    geoid.Geoid_Initialized = 1;


    """
    if(abs(coord_step) < 1.0e-10) coord_step = 99999.0
    if(abs(altitude_step_size) < 1.0e-10) altitude_step_size = 99999.0;
    if(fabs(time_step) < 1.0e-10) time_step = 99999.0;
    """

    cdef int NumTerms = ((magnetic_model.nMax + 1) * (magnetic_model.nMax + 2) / 2)
    cdef MAGtype_MagneticModel *TimedMagneticModel = MAG_AllocateModelMemory(NumTerms)
    cdef MAGtype_MagneticResults MagneticResultsSph, MagneticResultsGeo, MagneticResultsSphVar, MagneticResultsGeoVar
    cdef MAGtype_SphericalHarmonicVariables *SphVariables = MAG_AllocateSphVarMemory(magnetic_model.nMax)
    cdef MAGtype_GeoMagneticElements GeoMagneticElements
    cdef MAGtype_LegendreFunction *LegendreFunction = MAG_AllocateLegendreFunctionMemory(NumTerms)

    cdef MAGtype_CoordGeodetic coord
    cdef MAGtype_CoordSpherical coord_spherical
    cdef MAGtype_Date date

    axis_x = axis(minx, maxx, dx)
    axis_y = axis(miny, maxy, dy)
    axis_z = axis(minz, maxz, dz)
    axis_t = axis(mint, maxt, dt)

    cdef np.ndarray result = np.empty((
        axis_z.size, axis_y.size, axis_x.size, axis_t.size
        ), dtype=np.float64
    )

    for z in axis_z:
        for y in axis_y:
            for x in axis_x:
                coord.lambda_ = x
                coord.phi = y
                coord.HeightAboveGeoid = z

                if geoid.UseGeoid:
                    MAG_ConvertGeoidToEllipsoidHeight(&coord, &geoid)
                else:
                    coord.HeightAboveEllipsoid = coord.HeightAboveGeoid

                MAG_GeodeticToSpherical(ellipsoid, coord, &coord_spherical);
                MAG_ComputeSphericalHarmonicVariables(ellipsoid, coord_spherical, magnetic_model.nMax, SphVariables)
                MAG_AssociatedLegendreFunction(coord_spherical, magnetic_model.nMax, LegendreFunction)

                for d in axis_t:
                    date.DecimalYear = d

                    MAG_TimelyModifyMagneticModel(date, magnetic_model, TimedMagneticModel)
                    MAG_Summation(LegendreFunction, TimedMagneticModel, SphVariables[0], coord_spherical, &MagneticResultsSph)
                    MAG_SecVarSummation(LegendreFunction, TimedMagneticModel, SphVariables[0], coord_spherical, &MagneticResultsSphVar)
                    MAG_RotateMagneticVector(coord_spherical, coord, MagneticResultsSph, &MagneticResultsGeo)
                    MAG_RotateMagneticVector(coord_spherical, coord, MagneticResultsSphVar, &MagneticResultsGeoVar)
                    MAG_CalculateGeoMagneticElements(&MagneticResultsGeo, &GeoMagneticElements)
                    MAG_CalculateSecularVariationElements(MagneticResultsGeoVar, &GeoMagneticElements)
                    
                    if product_id == 1:
                        value = GeoMagneticElements.Decl
                    elif product_id == 2:
                        value = GeoMagneticElements.Incl
                    elif product_id == 3:
                        value = GeoMagneticElements.F
                    elif product_id == 4:
                        value = GeoMagneticElements.H
                    elif product_id == 5:
                        value = GeoMagneticElements.X
                    elif product_id == 6:
                        value = GeoMagneticElements.Y
                    elif product_id == 7:
                        value = GeoMagneticElements.Z
                    elif product_id == 8:
                        value = GeoMagneticElements.GV
                    elif product_id == 9:
                        value = GeoMagneticElements.Decldot
                    elif product_id == 10:
                        value = GeoMagneticElements.Incldot
                    elif product_id == 11:
                        value = GeoMagneticElements.Fdot
                    elif product_id == 12:
                        value = GeoMagneticElements.Hdot
                    elif product_id == 13:
                        value = GeoMagneticElements.Xdot
                    elif product_id == 14:
                        value = GeoMagneticElements.Ydot
                    elif product_id == 15:
                        value = GeoMagneticElements.Zdot
                    elif product_id == 16:
                        value = GeoMagneticElements.GVdot

                    result[
                        axis_z.index(z), axis_y.index(y), 
                        axis_x.index(x), axis_t.index(d)
                    ] = value
                    

    MAG_FreeMagneticModelMemory(TimedMagneticModel)
    MAG_FreeLegendreMemory(LegendreFunction)
    MAG_FreeSphVarMemory(SphVariables)

    return result


cdef class axis:
    cdef double minv
    cdef double maxv
    cdef double dv
    def __cinit__(self, double minv, double maxv, double dv):
        self.minv = minv
        self.maxv = maxv
        self.dv = dv

    @property
    def size(self):
        return int((self.maxv - self.minv) / self.dv) + 1

    def index(self, value):
        return int((value - self.minv) / self.dv)

    def __iter__(self):
        cdef double r = self.minv
        while r <= self.maxv:
            yield r
            r += self.dv


def to_year_fraction(v):
    def sinceEpoch(v): # returns seconds since epoch
        return time.mktime(v.timetuple())
    s = sinceEpoch

    year = v.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(v) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return v.year + fraction
