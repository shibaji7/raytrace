/*
 ==============================================================
                ground_fs_loss_matlab_wrapper.for
 ==============================================================
  Matlab wrapper to forward_scatter_loss.for which calculates the 
  the ground forward scatter power loss (dB) of the radio-waves.
       
  This is a MEX-file for MATLAB.
      
  Change history:
  11/09/06  M.A.Cervera    V1.0  Author
  10/05/07  M.A.Cervera    V2.0  Converted from Fortran
  18/02/10  M.A. Cervera   V2.1  Improved error handling
  13/07/17  D.J. Netherway V2.2  Allow array inputs
 
  13/02/22  M. A. Cervera  V2.2 
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

 ==============================================================
*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

int check_ref_data(char *files_to_check);

double forward_scatter_loss_(double *lat, double *lon, double *elev, double *freq);

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int  ii, mrows, ncols, arrsize, out_arrsize, mrows0, ncols0, ray_idx;
  char land_type;
  char error_mesg[200], filename[200], *refdata_dir;
  FILE *fopen(), *file_ptr;

  /* Input parameters to land_type.for */
  double lat, lon;    /* Geographic lat/lon of start-point */
  double elev;        /* forwared scattered elevation of ray at ground */
  double freq;        /* radio frequency of ray */
  
  /* Output parameters from forward_scatter_loss.for */
  double fs_loss;

  /* Other declarations */
  double origin_lat, origin_long, range, azim;

  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    mexErrMsgTxt("ERROR: 4 inputs required when calling ground_fs_loss.");
    return;
  }
  if (nlhs != 1) {
    mexErrMsgTxt("ERROR: 1 output required when calling ground_fs_loss.");
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
    arrsize = mrows*ncols;
    if (ii == 0) {
      mrows0 = mrows;
      ncols0 = ncols;
    }
    if (mrows != mrows0 || ncols != ncols0) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect - it should match first array size.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* Check to see if the reference data files exist. */
  if (check_ref_data("land_sea") == 0) return;

  /* Create matricies for the return arguments. */
  out_arrsize = mrows0 * ncols0;
  plhs[0] = mxCreateDoubleMatrix(mrows0, ncols0, mxREAL);

  /* loop over rays and call the computational routine */
  for (ray_idx = 0; ray_idx < out_arrsize; ray_idx++) {

    /* Copy the RHS pointers (input from MATLAB) to the input variables of
       forward_scatter_loss_ */
    lat  = *(mxGetPr(prhs[0]) + ray_idx);      
    lon  = *(mxGetPr(prhs[1]) + ray_idx);
    elev = *(mxGetPr(prhs[2]) + ray_idx);
    freq = *(mxGetPr(prhs[3]) + ray_idx);

    /*  Call the computational subroutine. */
    fs_loss = forward_scatter_loss_(&lat, &lon, &elev, &freq);

    /* Load the output data from firic into the LHS pointers, which are the 
       outputs to MATLAB. */
     *(mxGetPr(plhs[0])+ray_idx) = fs_loss;
  }

}
