#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


#include "GeomagnetismHeader.h"
#include "EGM9615.h"
//#include "GeomagnetismLibrary.c"

/*


WMM grid program.

The Geomagnetism Library is used to make a command prompt program. The program prompts
the user to enter a location, performs the computations and prints the results to the
standard output. The program expects the files GeomagnetismLibrary.c, GeomagnetismHeader.h,
WMM2010SHDF.COF and EGM9615.h to be in the same directory. 

Manoj.C.Nair@Noaa.Gov
April 21, 2011
 *  Revision Number: $Revision: 843 $
 *  Last changed by: $Author: awoods $
 *  Last changed on: $Date: 2012-04-25 16:23:56 -0600 (Wed, 25 Apr 2012) $

 */



int main()
{
    MAGtype_MagneticModel * MagneticModels[1];
    MAGtype_Ellipsoid Ellip;
    MAGtype_CoordGeodetic minimum, maximum;
    MAGtype_Geoid Geoid;
    MAGtype_Date startdate, enddate;
    int ElementOption, PrintOption, i, epochs = 1;
    double cord_step_size, altitude_step_size, time_step_size;
    char filename[] = "WMM.COF";
    char OutputFilename[32];
    char VersionDate_Large[] = "$Date: 2012-04-25 16:23:56 -0600 (Wed, 25 Apr 2012) $";
    char VersionDate[12];

    MAG_robustReadMagModels(filename, &MagneticModels, 1);

    strncpy(VersionDate, VersionDate_Large + 39, 11);
    VersionDate[11] = '\0';

    MAG_SetDefaults(&Ellip, &Geoid);
    //MAG_InitializeGeoid(&Geoid);
    /* Set EGM96 Geoid parameters */
    Geoid.GeoidHeightBuffer = GeoidHeightBuffer;
    Geoid.Geoid_Initialized = 1;
    /* Set EGM96 Geoid parameters END */
    //MAG_GeomagIntroduction_WMM(MagneticModels[0]);  /* Print out the WMM introduction */
    printf("\n\n Welcome to the World Magnetic Model (WMM) C-Program\nof the US National Geophysical Data Center\n\t\t--- Grid Calculation Program ----\n\t       --- Model Release Date: 15 Dec 2009 ---\n\t     --- Software Release Date: %s ---\n", VersionDate);
    printf("This program may be used to generate a grid/volume of magnetic field values\n");
    printf("over latitude, longitude, altitude and time axes. To skip an axis,\n");
    printf("keep the start and end values the same and enter zero for the step size.\n");

    printf("\n\n                     Enter grid parameters \n\n");


    /* Get the Lat/Long, Altitude, Time limits from a user interface and print the grid to screen */

    MAG_GetUserGrid(&minimum, &maximum, &cord_step_size, &altitude_step_size, &time_step_size, &startdate, &enddate, &ElementOption, &PrintOption, OutputFilename, &Geoid);
    MAG_Grid(minimum, maximum, cord_step_size, altitude_step_size, time_step_size, MagneticModels[0], &Geoid, Ellip, startdate, enddate, ElementOption, PrintOption, OutputFilename);

    for(i = 0; i < epochs; i++) MAG_FreeMagneticModelMemory(MagneticModels[i]);




    printf("\nPress any key to exit...\n");
    getchar();

    return 0;
}
