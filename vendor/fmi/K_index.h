/****************************************************************************/
/*                                                                          */
/*                                K_index.h                                 */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/* This is a C language header file that defines the routines for computing */
/* daily K indices by the FMI-method.                                       */
/*                                                                          */
/* There is one routine that the user needs to call in order to compute     */
/* K indices. That function is called 'K_index' and its definition is :     */
/*                                                                          */
/* void K_index(long *X_data,long *Y_data,long SampleStep,long K9_limit,    */
/*              long Longitude,long MissingData,long *K_table)              */
/*                                                                          */
/* The parameters are:                                                      */
/*      X_data      Pointer to the start of the data array containing the   */
/*                  X-component of the field values. The field values are   */
/*                  stored as long's and may be whatever unit one uses      */
/*                  (e.g. 1 nT or 0.1 nT). If one normally uses decimal     */
/*                  points (e.g. 49123.5) the values must be converted to   */
/*                  long's (e.g. 491235) before the routine can be called.  */
/*                  In order to compute the K indices for one day the FMI   */
/*                  method needs data for previous day and next day.        */
/*                  Therefore the X_data array must contain data exactly    */
/*                  for three days (i.e if sampling rate is one minute the  */
/*                  X_data array must have 3*1440 = 4320 elements).         */
/*      Y_data      Pointer to the start of the data array containing the   */
/*                  Y-component of the field values.                        */
/*      SampleStep  Interval between successive data points in seconds      */
/*                  (e.g. if minute values then SampleStep = 60).           */
/*      K9_limit    K=9 limit for the particular observatory. Here the      */
/*                  unit is the same as used in the field values. (e.g if   */
/*                  the K=9limit is 750 nT for the observatory and field    */
/*                  values are expressed with 0.1 nT accuracy then          */
/*                  K9_limit = 7500).                                       */
/*      Longitude   This is the integral part of the longitude of the       */
/*                  observatory whose K indices are to be computed. The FMI */
/*                  method uses it to determine the time of local midnight. */
/*                  (e.g. Nurmijarvi longitude is 24.65, so Longitude=25).  */
/*      MissingData Marker for missing data point (e.g. 999999).            */
/*      K_table     This is the table where the computed eight K values are */
/*                  returned. So one should define in his/her program an    */
/*                  array: long K[8]. If the program is unable to compute a */
/*                  K value (due to missing data) a value -1 is substituted */
/*                  for that particular element.                            */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/*                            Lasse Hakkinen                                */
/*                    Finnish Meteorological Institute                      */
/*                        Department of Geophysics                          */
/*                              P.O.Box 503                                 */
/*                      FIN-00101, Helsinki, Finland                        */
/*                      e-mail: Lasse.Hakkinen@fmi.fi                       */
/*                      phone : (+358)-9-19294634                           */
/*                      fax   : (+358)-9-19294603                           */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/
/*  Version history:                                                        */
/*                                                                          */
/*  1.3  11.03.1999  Fixed a bug in ComputeMean routine where a potential   */
/*                   "long overflow" error might be produced if one uses    */
/*                   1 second data and uses 0.1 nT resolution.              */
/*                   This bug was pointed out by Patrik Johansson.          */
/*  1.11 20.10.1998  Slight cosmetic changes.                               */
/*  1.1  02.12.1992  Fixed a bug which caused erronouos results if the day  */
/*                   started with missing data or if the previous day ended */
/*                   with missing data. The bug was in the routine          */
/*                   'InterpolateMean'.                                     */
/*  1.0  29.09.1992  First C language version                               */
/****************************************************************************/

#include <stdlib.h>
#include <math.h>

/*--------------------------------------------------------------------------*/
/* Parameters of the FMI-method. These default values are best for several  */
/* observatories. Changing the values may in some cases result in slightly  */
/* better agreement between hand-scaled and computed K-indices.             */
/*--------------------------------------------------------------------------*/

#define NightLength  90*60  /* Night fitting length in seconds              */
#define DawnLength   60*60  /* Dawn fitting lengths in seconds              */
#define DuskLength   60*60  /* Dusk fitting lengths in seconds              */
#define HarmOrder    5      /* Order of harmonic fitting                    */
#define KExponent    3.3    /* Exponent in K^KExponent                      */


/*--------------------------------------------------------------------------*/
/* Other parameters used in the analysis.                                   */
/*--------------------------------------------------------------------------*/

#define MaxGapLength 15*60 /* Maximum length of gap of missing data points  */
                           /* in seconds. Smaller gaps will be interpolated */
                           /* Larger gaps will lead to missing K-values.    */
#define PI 3.141592654


/*--------------------------------------------------------------------------*/
/* The global variables used in the computation of the K indices.           */
/*--------------------------------------------------------------------------*/

long PointsInDay;       /* Number of datapoints in one day                  */
long DataCount;         /* Number of datapoints in three days               */
long StepSize;          /* Data sampling interval in seconds                */
long MissingPoint;      /* Marker for missing data point                    */
long *X_temp;           /* Pointer to a temporary data array for X-comp.    */
long *Y_temp;           /* Pointer to a temporary data array for Y-comp.    */
long *X_harm;           /* Pointer to the harmonic fitting array of X-curve */
long *Y_harm;           /* Pointer to the harmonic fitting array of Y-curve */
long K_limit[9];        /* Array of limiting values of the K indices        */
long DefLength[24];     /* Default fitting length for each hour             */
long K_length[10];      /* Fitting lengths for each K value                 */
long FitLength[24];     /* Computed fitting length for each hour            */
long X_mean[24];        /* Mean value for each hour over the fitting area   */
long Y_mean[24];        /* Mean value for each hour over the fitting area   */


/*--------------------------------------------------------------------------*/
/* Initialize the global variables declared above. Since a priori we don't  */
/* know the sampling rate and the size of the data arrays we have to        */
/* allocate space dynamically. This is done by using the malloc-function.   */
/*--------------------------------------------------------------------------*/

static void InitGlobals(long Step,long K9_limit,long Longitude,long MissingData)
{
    long i;              /* Dummy index                                     */
    long BlockSize;      /* Size (in bytes) of the three day data           */
    long DummyTable[24]; /* A table used in adjusting the DefLength table   */
    long DefLimit[9] = {75,150,300,600,1050,1800,3000,4950,7500};
                         /* Limits for different K indices (unit = 0.1 nT)  */

    PointsInDay  = (24*3600)/Step;
    DataCount    = 3*PointsInDay;
    StepSize     = Step;
    MissingPoint = MissingData;
    BlockSize    = sizeof(long)*PointsInDay;
    X_temp = (long *) malloc(BlockSize);    /* Allocate space for tables */
    Y_temp = (long *) malloc(BlockSize);
    X_harm = (long *) malloc(BlockSize);
    Y_harm = (long *) malloc(BlockSize);

    /* Set the limit values for K indices */
    for (i=0;i<9;i++) K_limit[i] = (K9_limit*DefLimit[i])/7500;

    /* Set the default fitting lenghts for each hour */
    for (i=0;i<3;i++)    DefLength[i] = NightLength;
    for (i=3;i<6;i++)    DefLength[i] = DawnLength;
    for (i=6;i<18;i++)   DefLength[i] = 0;
    for (i=18;i<21;i++)  DefLength[i] = DuskLength;
    for (i=21;i<24;i++)  DefLength[i] = NightLength;

    /* Rotate the DefLength table to adjust the midnight */
    for (i=0;i<24;i++) DummyTable[i] = DefLength[(i+Longitude/15) % 24];
    for (i=0;i<24;i++) DefLength[i] = DummyTable[i];

    /* Compute the fitting lengths for different K-values */
    K_length[0] = 0;
    for (i=1;i<10;i++) {
        K_length[i] = (long) (60*exp(KExponent*log(i)));
        if (K_length[i] > 18*3600) K_length[i] = 18*3600;/* max = 18 hours */
    }
}


/*--------------------------------------------------------------------------*/
/* Go through the three-day data and interpolate gaps whose lengths are     */
/* smaller than MaxGapLength. If gaps are larger than this leave them       */
/* untouched.  K index computation routine will notice the gaps and return  */
/* a Missing K value for the particular three hour period.                  */
/*--------------------------------------------------------------------------*/

static void FillGaps(long *Data_array)
{
    long  i,j;              /* Indices for data points                  */
    long  GapIndex;         /* Index of the first missing data point    */
    long  GapLength;        /* Number of missing data points            */
    long  Value1;           /* Last non-missing data point before gap   */
    long  Value2;           /* First non-missing data point after gap   */

    i = 0;
    while (i < DataCount) {            /* Go through all the three days */
        /** Find the next missing value **/
        while ((i < DataCount) && (Data_array[i] != MissingPoint)) i++;

        /** If a gap was found handle it **/
        if (i<DataCount) {
            GapIndex = i;
            /* Find the last missing point */
            while ((i < DataCount) && (Data_array[i] == MissingPoint)) i++;

            /** Interpolate the gap if possible. We cannot extrapolate **/
            if ((GapIndex > 0) && (i < DataCount)) {
                /** Now i is the index of first non-missing value    **/
                /** and GapIndex is the index of first missing value **/
                GapLength = i-GapIndex;
                /** Interpolate if the gap is small enough **/
                if (GapLength*StepSize < MaxGapLength) {
                    Value1 = Data_array[GapIndex-1];
                    Value2 = Data_array[i];
                    for (j=1;j<=GapLength;j++)
                        Data_array[GapIndex++] =
                                    Value1+(j*(Value2-Value1))/(GapLength+1);
                }
            }
        }
    }
}

/*--------------------------------------------------------------------------*/
/* Copy one day data from FromArray to ToArray                              */
/*--------------------------------------------------------------------------*/

static CopyData(long *FromArray, long *ToArray)
{
    long i;

    FromArray += PointsInDay;   /* FromArray contains data for three days.  */
                                /* Select the middle one.                   */
    for (i=0;i<PointsInDay;i++) *ToArray++ = *FromArray++;
}

/*--------------------------------------------------------------------------*/
/* Compute the difference : DiffArray = Array1 - Array2.                    */
/*--------------------------------------------------------------------------*/

static ComputeDiff(long *Array1, long *Array2, long *DiffArray)
{
    long i;

    Array1 += PointsInDay;      /* Array1 contains data for three days.     */
                                /* Select the middle one.                   */
    for (i=0;i<PointsInDay;i++) {
        if ((*Array1 != MissingPoint) && (*Array2 != MissingPoint))
            *DiffArray++ = *Array1++ - *Array2++;
        else {
            *DiffArray++ = MissingPoint;
            Array1++;
            Array2++;
        }
    }
}

/*--------------------------------------------------------------------------*/
/* Compute the difference between the maximum and minimum value in a three  */
/* hour part in the Data_array and determine the corresponding K index.     */
/*--------------------------------------------------------------------------*/

static long K_MaxMin(long *Data_array, long HourIndex)
{
    long i;             /* Dummy index                  */
    long Max,Min;       /* Maximum and minimum values   */

    /* Find the maximum and minimum values */
    Data_array += HourIndex*PointsInDay/8;  /* Find the right 3-hour block */
    Max = *Data_array;
    Min = *Data_array;
    for (i=1;i<PointsInDay/8;i++) {
        if (*Data_array > Max) Max = *Data_array; else
        if (*Data_array < Min) Min = *Data_array;
        Data_array++;
    }

    /* Find the K index corresponding to Max-Min */
    if (Max == MissingPoint) return (-1);   /* -1 if unable to compute K */
    else {      /* Go through the KLimit-table */
        i = 0;
        while ((i<9) && (Max-Min > K_limit[i])) i++;
        return (i);
    }
}


/*--------------------------------------------------------------------------*/
/* Determine K indices by using the max-min method either directly to the   */
/* original data (Subtract = false) or to the difference between original   */
/* data and harmonic fitting (Subtract = true). The final K index is the    */
/* larger one from the computed K's for the X and Y components.             */
/*--------------------------------------------------------------------------*/

static void FillKTable(long *K_table,long *X_data,long *Y_data,long Subtract)
{
    long i;             /* Index variable                           */
    long K_X[8];        /* Three hour K indices from X-component    */
    long K_Y[8];        /* Three hour K indices from Y-component    */

    /* Fill first the X_temp and Y_temp arrays */
    if (Subtract) {
        ComputeDiff(X_data,X_harm,X_temp);
        ComputeDiff(Y_data,Y_harm,Y_temp);
    }
    else {
        CopyData(X_data,X_temp);
        CopyData(Y_data,Y_temp);
    }

    /* Now compute K indices for X and Y from max-min values */
    /* for each three-hour interval and take the larger one  */
    /* as the final K index.                                 */
    for (i=0;i<8;i++) {
        K_X[i] = K_MaxMin(X_temp,i);
        K_Y[i] = K_MaxMin(Y_temp,i);
        if (K_X[i] > K_Y[i]) K_table[i] = K_X[i]; else K_table[i] = K_Y[i];
    }
}


/*--------------------------------------------------------------------------*/
/* Find the fitting lengths for each hour.                                  */
/*--------------------------------------------------------------------------*/

static void FindLengths(long *K_table)
{
    long i;

    for (i=0;i<24;i++) {                /* go through all hours */
        if (K_table[i/3] < 0) FitLength[i] = 30*60+DefLength[i];
            else FitLength[i] = 30*60+DefLength[i]+K_length[K_table[i/3]];
        if (FitLength[i] > 24*3600)     /* not more than one day */
            FitLength[i] = 24*3600;
    }
}

/*--------------------------------------------------------------------------*/
/* Compute mean values for each hour over the fitting interval.             */
/*--------------------------------------------------------------------------*/

static long ComputeMean(long hour, long *Data_array)
{
    long i;           /* Dummy index                                        */
    long MiddlePoint; /* Index of the middle point of the fitting interval  */
    long WingLength;  /* Number of datapoints before and after middle point */
    long Count;       /* Number of datapoints in the fitting interval       */
    double Sum = 0.0; /* Sum of field values                                */

    MiddlePoint = PointsInDay+(60*(60*hour+30))/StepSize;
    WingLength  = FitLength[hour]/StepSize;
    Data_array += MiddlePoint-WingLength;
    Count = 2*WingLength+1;

    for (i=0;i<Count;i++) {
        if (*Data_array == MissingPoint) return (MissingPoint);
        else Sum += *Data_array++;
    }
    return (long)(Sum/Count);
}

/*--------------------------------------------------------------------------*/
/* Interpolate the missing data values in the hourly mean value tables.     */
/*--------------------------------------------------------------------------*/

static void InterpolateMean(long *Mean_array)
{
    long i,hour;
    long hour0 = 0;
    long hour1 = 23;
    long GapLength;         /* Number of missing data points                */
    long Value1;            /* Last non-missing data point before gap       */
    long Value2;            /* First non-missing data point after gap       */
    long found;             /* Boolean flag indicating missing data points  */

    /* Extrapolate the first missing points */
    while ((hour0 < 24) && (Mean_array[hour0] == MissingPoint)) hour0++;
    if (hour0 < 24)
        for (hour=0;hour<hour0;hour++) Mean_array[hour] = Mean_array[hour0];

    /* Extrapolate the last missing points */
    while ((hour1 > hour0) && (Mean_array[hour1] == MissingPoint)) hour1--;
    for (hour = hour1+1;hour<24;hour++) Mean_array[hour] = Mean_array[hour1];

    /* Interpolate the missing points between hour0 and hour1 */
    /* Both hour0 and hour1 are hours where data exists       */
    while (hour0 < hour1) {
        do {
            found = (Mean_array[++hour0] == MissingPoint);
        } while ((hour0 < hour1) && (!found));
        if (found) {
            hour = hour0;           /* first missing hour */
            while ((hour0<hour1) && (Mean_array[hour0]==MissingPoint))
                hour0++;
            GapLength = hour0-hour;
            Value1 = Mean_array[hour-1];
            Value2 = Mean_array[hour0];
            for (i=1;i<=GapLength;i++)
                Mean_array[hour++] = Value1+(i*(Value2-Value1))/(GapLength+1);
        }
    }
}


/*--------------------------------------------------------------------------*/
/* Compute the Hour averages and interpolate missing hour values if needed. */
/*--------------------------------------------------------------------------*/

static void ComputeHourAverages(long *X_data, long *Y_data)
{
    long hour;

    for (hour=0;hour<24;hour++) {
        X_mean[hour] = ComputeMean(hour,X_data);
        Y_mean[hour] = ComputeMean(hour,Y_data);
    }
    
    /* Interpolate or extrapolate the possible gaps */
    InterpolateMean(X_mean);
    InterpolateMean(Y_mean);
}




/*--------------------------------------------------------------------------*/
/* Compute the harmonic fit to the Hour-averages and put data values into   */
/* X_harm- and Y_harm-tables.                                               */
/*--------------------------------------------------------------------------*/

static void ComputeHarmonicFit()
{
    long i,k;           /* Dummy indices */
    float X_coeff;
    float Y_coeff;
    long *Xptr,*Yptr;
    float ReX[13],ImX[13];
    float ReY[13],ImY[13];
    float si,co;
    float angle,angle2;
    long  t,t0,t1;

    t0 = 30*60;         /* Middle point of first hour in seconds */
    t1 = 23*3600+t0;    /* Middle point of last hour in seconds  */

    /* Rotate the curve so that start and end points are equal */
    X_coeff = ((float)(X_mean[23]-X_mean[0]))/(t1-t0);
    Y_coeff = ((float)(Y_mean[23]-Y_mean[0]))/(t1-t0);
    t = t0;
    for (i=0;i<24;i++) {
        X_mean[i] -= (long)(X_coeff*(t-t0));
        Y_mean[i] -= (long)(Y_coeff*(t-t0));
        t += 3600;
    }

    /* Compute the Fourier coefficients */
    for (i=0;i<=HarmOrder;i++) {
        ReX[i] = X_mean[0];
        ImX[i] = 0.0;
        ReY[i] = Y_mean[0];
        ImY[i] = 0.0;
        angle = -i*(2.0*PI/24.0);
        for (k=1;k<24;k++) {
            si = sin(k*angle);
            co = cos(k*angle);
            ReX[i] += X_mean[k]*co;
            ImX[i] += X_mean[k]*si;
            ReY[i] += Y_mean[k]*co;
            ImY[i] += Y_mean[k]*si;
        }
    }

    /* Compute the inverse fourier transform taking into account */
    /* only terms up to HarmOrder.                               */
    Xptr = X_harm;
    Yptr = Y_harm;
    t = 0;
    angle = 2.0*PI*(23.0/24.0)/(t1-t0);
    for (i=0;i<PointsInDay;i++) {
        *Xptr = (long)ReX[0];
        *Yptr = (long)ReY[0];
        angle2 = (t-t0)*angle;
        for (k=1;k<=HarmOrder;k++) {
            si = sin(k*angle2);
            co = cos(k*angle2);
            *Xptr += (long)(2.0*(ReX[k]*co-ImX[k]*si));
            *Yptr += (long)(2.0*(ReY[k]*co-ImY[k]*si));
        }
        /* Rotate the curve back */
        *Xptr = (long)(*Xptr/24.0+X_coeff*(t-t0));
        *Yptr = (long)(*Yptr/24.0+Y_coeff*(t-t0));
        Xptr++;
        Yptr++;
        t += StepSize;
    }
}


/*--------------------------------------------------------------------------*/
/* Compute the Ak number for one day using the already computed K indices.  */
/*--------------------------------------------------------------------------*/

long ComputeAk(long *K_table)
{
    static long AkTable[10] = {0,3,7,15,27,48,80,140,240,400};
    long Ak     = 0;            /* Weighted sum of K indices */
    long count  = 0;            /* Number of valid K indices */
    long i;

    for (i=0;i<8;i++) {
        if (K_table[i] >= 0) {
            Ak += AkTable[K_table[i]];
            count++;
        }
    }
    if (count > 0) return ((Ak+count/2)/count); /* Round to nearest integer */
    else return (-1);
}

/*--------------------------------------------------------------------------*/
/* Release the space allocated to temporary arrays.                         */
/*--------------------------------------------------------------------------*/

static void FreeArrays()
{
    free(X_temp);
    free(Y_temp);
    free(X_harm);
    free(Y_harm);
}

/*--------------------------------------------------------------------------*/
/* Compute the K indices for one day. This is the only function that the    */
/* user should ever need to call in this header file.                       */
/*--------------------------------------------------------------------------*/

void K_index(long *X_data,long *Y_data,long SampleStep,
             long K9_limit,long Longitude,long MissingData,long *K_table)
{
    /* Initialize global variables and allocate space to temporary arrays   */
    InitGlobals(SampleStep,K9_limit,Longitude,MissingData);

    FillGaps(X_data);   /* Try to interpolate possible gaps in the data */
    FillGaps(Y_data);

    /*** Step 1: Compute preliminary K indices by max-min method ***/
    FillKTable(K_table,X_data,Y_data,0); /* Here 0 implies not to subtract  */
                                         /* harmonic fitting from orig data */

    /*** Step 2: Compute the fitting lengths and find average values ***/
    /*** for each hour                                               ***/
    FindLengths(K_table);
    ComputeHourAverages(X_data,Y_data);

    /*** Step 3: Compute the harmonic fit to hour averages ***/
    ComputeHarmonicFit();

    /*** Step 4: Compute new K values from Original data - Harmonic fit ***/
    FillKTable(K_table,X_data,Y_data,1); /* Here 1 implies to subtract      */
                                         /* harmonic fitting from orig data */

    /*** Step 5: Compute the fitting lengths and find average values ***/
    /*** for each hour                                               ***/
    FindLengths(K_table);
    ComputeHourAverages(X_data,Y_data);

    /*** Step 6: Compute the harmonic fit to hour averages ***/
    ComputeHarmonicFit();

    /*** Step 7: Compute final K values from Original data - Harmonic fit ***/
    FillKTable(K_table,X_data,Y_data,1);

    FreeArrays();       /* Release the space allocated to temporary arrays */
}

