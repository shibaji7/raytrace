/*
 ================================================================
                   nrlmsise00_matlab_wrapper.c
 ================================================================
  Matlab wrapper for the NRLMSISE00 subroutines distributed in iri2012
       
  This is a MEX-file for MATLAB.
      
  Change history:
  30/08/12  V1.0  L.H.Pederick
     Author

  07/09/2018 V1.1  M.A. Cervera 
     Now accepts F10.7 and Ap index as inputs. R12 is no longer used.
     If F10.7 and Ap are no supplied then these values are read using
     the routine apfmsis supplied by IRI2016.

  26/07/2019 V1.2 M.A. Cervera 
     Fixed bug: The local apparent solar time, stl, was being rounded down 
     to the nearest hour due to a type casting error. This had the effect of 
     nrlmsise00 returning its output on the hour instead of hour and minute. 

  13/02/2022  M. A. Cervera  V1.3 
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

 ================================================================
*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

void nrlmsise00_calc_(int *iyd, float *sec, float *alt, float *glat, 
           float *glong, float *stl, float *F107_81, float *F107_prior_day,
	   float *Ap_daily, int *mass, float *d, float *t);

void apfmsis_call_(int *UT, float *F107_day, float *F107_prior_day, 
                   float *F107_81, float *F107_365, float *Ap);

int  julday_(int *day, int* month, int *year);

int check_ref_data(char *files_to_check);

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  float sec, alt, glat, glong, stl, Ap_daily;
  float F107_day, F107_prior_day, F107_81, F107_365, Ap[7], d[9], t[2];
  int UT[5], iyd, mass, nmonth;
  int year, month, day, hour, min, days_in_year, day_of_year, dd1, mm1, dd2,
      mm2;
  int ii, jj, mrows, ncols, arrsize;
  char error_mesg[200];


  /* Check for proper number of arguments. */
  if ((nrhs != 4) && (nrhs != 7)) {
    mexErrMsgTxt("ERROR: either 4 or 7 inputs required when calling nrlmsise00.");
    return;
  }

  if (nlhs != 2) {
    mexErrMsgTxt("ERROR: 2 outputs required when calling nrlmsise00.");
    return;
  }

  /* Check to ensure the inputs are numeric (not strings). */
  for (ii = 0; ii < nrhs; ii++) {
    if (mxIsNumeric(prhs[ii]) == 0) {
      snprintf(error_mesg, 200, "ERROR: Input %d must be numeric.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* Get and check the size of the input arrays. */
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if (mrows*ncols != 5) {
    mexErrMsgTxt("ERROR: Input 4 array size is incorrect (must be 5x1).");
    return;
  }

  if (nrhs >= 5) {
    for (ii = 4; ii < 7; ii++) { 
      mrows = mxGetM(prhs[ii]);
      ncols = mxGetN(prhs[ii]);
      if (mrows*ncols != 1) {
        snprintf(error_mesg, 200, "ERROR: Input %d must be a scalar.", ii+1);
        mexErrMsgTxt(error_mesg);
        return;
      }
    }
  }

  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  arrsize = mrows*ncols;

  for (ii = 0; ii < 3; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    if (mrows*ncols != arrsize) {
      mexErrMsgTxt("ERROR: Inputs 1 to 3 must have the same array size.");
    }
  }

  /* Check to see if the reference data files exist. */
  if (check_ref_data("iri2016") == 0) return;

  /* Copy the RHS pointers (input from MATLAB) to the input variables of
     ther igrf routines */
  for (ii = 0; ii < 5; ii++) {
    UT[ii] = *(mxGetPr(prhs[3]) + ii);
  }

  /* Calculate inputs to nrlmsise00_calc */
  year = UT[0];
  month = UT[1];
  day = UT[2];
  hour = UT[3];
  min = UT[4];
  
  dd1 = 0; mm1 = 1; dd2 = 31; mm2 = 12;
  day_of_year  = julday_(&day, &month, &year) - julday_(&dd1, &mm1, &year);
  days_in_year = julday_(&dd2, &mm2, &year) - julday_(&dd1, &mm1, &year);
  
  iyd = (year % 100)*1000 + day_of_year;
  sec = min*60 + hour*60*60;

  mass = 48; /* calculate all gases */

 
  /* Get the F10.7 parameters and Ap magnetic index. If not supplied */
  /* then read in F10.7 parameters and Ap for the supplied UT using apfmsis */
  if (nrhs == 7) {
    F107_prior_day = mxGetScalar(prhs[4]);
    F107_81 = mxGetScalar(prhs[5]);
    Ap_daily = mxGetScalar(prhs[6]);
  } else {
    apfmsis_call_(UT, &F107_day, &F107_prior_day, &F107_81, &F107_365, Ap);
    Ap_daily = Ap[1];  
  }

  /* Setup output arrays */
  plhs[0] = mxCreateNumericMatrix(9, arrsize, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(2, arrsize, mxDOUBLE_CLASS, mxREAL);

  if (arrsize == 1) {
    glat = mxGetScalar(prhs[0]);
    glong = mxGetScalar(prhs[1]);
    alt = mxGetScalar(prhs[2]);
    stl = ((float) hour + ((float) min)/60 + glong/15); /* local solar time */
    if (stl > 24.0) stl -= 24.0;
    if (stl < 0.0) stl += 24.0;

    nrlmsise00_calc_(&iyd, &sec, &alt, &glat, &glong, &stl, &F107_81, 
                     &F107_prior_day, &Ap_daily, &mass, d, t);
    for( ii = 0; ii < 9; ii++ )
        (mxGetPr(plhs[0]))[ii] = d[ii];

    for( ii = 0; ii < 2; ii++ )
        (mxGetPr(plhs[1]))[ii] = t[ii];
  } else {
      for ( jj = 0; jj < arrsize; jj++ ) {
        glat = *(mxGetPr(prhs[0]) + jj);
        glong = *(mxGetPr(prhs[1]) + jj);
        alt = *(mxGetPr(prhs[2]) + jj);
        stl = ((float) hour + ((float) min)/60 + glong/15); /* local solar time */
        if (stl > 24.0) stl -= 24.0;
        if (stl < 0.0) stl += 24.0;
        nrlmsise00_calc_(&iyd, &sec, &alt, &glat, &glong, &stl, &F107_81, 
                         &F107_prior_day, &Ap_daily, &mass, d, t);
        
        for( ii = 0; ii < 9; ii++ )
            (mxGetPr(plhs[0]))[ii + jj*9] = d[ii];

        for( ii = 0; ii < 2; ii++ )
            (mxGetPr(plhs[1]))[ii + jj*2] = t[ii];
    }
  }
}
