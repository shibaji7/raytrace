/*
 ==============================================================
                     dop_spread_eq_matlab_wrapper.c
 ==============================================================
  Matlab wrapper for dop_spread_eq.for
      
  This creates a MEX-file for MATLAB.
     
  Change history:
  31/10/05  M.A.Cervera  V1.0  Author
  02/05/07  M.A.Cervera  V2.0  Converted from Fortran
  18/02/10  M.A. Cervera V2.1  Improved error handling

  13/02/22  M. A. Cervera  V2.2 
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

 ==============================================================
*/

#include "mex.h"
#include "../iono_structures.h"

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

int check_ref_data(char *files_to_check);

float dop_spread_eq_(float *glon, float *glat, float *year, int *month,
		     int *day, int *hour, int*minute, float *R12);

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int ii, mrows, ncols, arrsize, in_arrsize[5], out_arrsize;
  char error_mesg[200];

  /* Input parameters to cs.for */
  float glat, glon ;                 /* Geographic lat/long */
  int month, day, hour, minute;      /* UT */
  float year, R12, bearing;          /* R12 index and ray bearing */

  /* Output parameters from dop_spread_eq_ */
  float doppler_spread;

  /* Define the expected input array sizes */
  in_arrsize[0] = 1;
  in_arrsize[1] = 1;
  in_arrsize[2] = 5;
  in_arrsize[3] = 1;

  /* Define the expected output array sizes */
  out_arrsize = 1;

  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    mexErrMsgTxt("ERROR: 4 inputs required when calling dop_spread_eq.");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("ERROR: 1 output required when calling dop_spread_eq.");
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
  for (ii = 0; ii < nrhs; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* Check to see if the reference data files exist. */
  /* currently does not use reference data files so don't check */
  /* if (check_ref_data() == 0) return; */

  /* Copy the RHS pointers (input from MATLAB) to the input variables of cs_ */
  glat = (float) mxGetScalar(prhs[0]);
  glon = (float) mxGetScalar(prhs[1]);
  year = *(mxGetPr(prhs[2]) + 0);
  month = *(mxGetPr(prhs[2]) + 1);
  day = *(mxGetPr(prhs[2]) + 2);
  hour = *(mxGetPr(prhs[2]) + 3);
  minute = *(mxGetPr(prhs[2]) + 4);
  R12 = mxGetScalar(prhs[3]);

  /* Create matricies for the return arguments. */
  plhs[0] = mxCreateDoubleMatrix(1, out_arrsize, mxREAL);      

  /* Call the computational subroutine. */
  doppler_spread = dop_spread_eq_(&glon, &glat, &year, &month, &day, &hour,
				  &minute, &R12);

  /* Load the output data from cs into the LHS pointers, which are the 
  outputs to MATLAB.
  */
  *(mxGetPr(plhs[0])) = doppler_spread;

}
