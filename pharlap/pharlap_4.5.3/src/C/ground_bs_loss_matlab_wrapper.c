/*
  ==============================================================
                 ground_bs_loss_matlab_wrapper.c
  ==============================================================
   Matlab wrapper to land_type_ and subsequent calculation of 
   the power loss (dB) of the radio-waves back-scattered from the 
   ground.
       
   This is a MEX-file for MATLAB.
      
   Change history:
   11/09/06  M.A.Cervera  V1.0  Author
   10/05/07  M.A.Cervera  V2.0  Converted from Fortran
   18/02/10  M.A. Cervera V2.1  Improved error handling
   13/07/17  M.A. Cervera V2.2  Allow array inputs
 
   13/02/22  M. A. Cervera  V2.3 
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

  ==============================================================
*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

int check_ref_data(char *files_to_check);

void land_sea_(double* lat, double* lon, char *land_type); 

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int  ii, mrows, ncols, mrows0, ncols0, ray_idx;
  char land_type[2], error_mesg[200], filename[200], *refdata_dir;
  FILE *fopen(), *file_ptr;

  /* Input parameters to land_type_ */
  double lat,lon;               /* Geographic lat/lon of start-point */

  /* Output parameters from back_scatter_loss_ */
  double back_scatter_loss;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    mexErrMsgTxt("ERROR: 2 inputs required when calling ground_bs_loss.");
    return;
  }

  if (nlhs != 1) {
    mexErrMsgTxt("ERROR: 1 outputs required when calling ground_bs_loss.");
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
  for (ii = 0; ii < nrhs; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    if (ii == 0) {
      mrows0 = mrows;
      ncols0 = ncols;
    }
    if (mrows != mrows0 || ncols != ncols0) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect - it should match first array size", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* Check to see if the reference data files exist. */
  if (check_ref_data("land_sea") == 0) return;

  /*  Create matricies for the return arguments. */
  plhs[0] = mxCreateDoubleMatrix(mrows0, ncols0, mxREAL);

  /* loop over rays and call the computational routine */
  for (ray_idx = 0; ray_idx < mrows0*ncols0; ray_idx++) {

    /* Initialize the land_type string variable and call land_sea_ to determine 
       whether the point of interest is land or sea.*/
    land_type[0] = 'S';
    land_type[1] = '\0';
    back_scatter_loss = 0.0025;                            /* this is sea */

    /* Copy the RHS pointers (input from MATLAB) to the input variables of 
       land_sea_ */
    lat = *(mxGetPr(prhs[0]) + ray_idx);
    lon = *(mxGetPr(prhs[1]) + ray_idx);

    land_sea_(&lat, &lon, land_type);
    if (strcmp(land_type, "L") == 0) {                     /* this is land */
      back_scatter_loss = back_scatter_loss * 0.5; 
    }  
    back_scatter_loss = -10*log10(back_scatter_loss); 

    /* Load the output data from firic into the LHS pointers, which are the 
       outputs to MATLAB. */
    *(mxGetPr(plhs[0])+ray_idx) = back_scatter_loss;
  }

}
