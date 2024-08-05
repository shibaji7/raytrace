/*
 ==============================================================
                   abso_bg_matlab_wrapper.c
 ==============================================================
  Matlab wrapper to abso_bg.for which calculates the D-Region 
  absorption of the radio-waves using Bradley-George.
       
  This creats a MEX-file for MATLAB.
      
  Change history:
  11/09/06  M.A.Cervera    V1.0  Author
  20/04/07  M.A.Cervera    V2.0  Converted from Fortran
  11/09/09  M.A. Cervera   V2.1  Modified to reflect changes to abso_bg.for
  18/02/10  M.A. Cervera   V2.2  Improved error handling
  13/07/17  D.J. Netherway V2.3  Allow array inputs first four arguments
 
  13/02/22  M. A. Cervera  V2.4
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

double abso_bg_(double *lat, double *lon, double *elev, double *freq, 
                int *UT, double *R12_index, int *O_mode, int *warning_flag);

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, arrsize, out_arrsize, ii, arrsize0, mrows0, ncols0, ray_idx;
  int in_arrsize[8];
  char error_mesg[200];

  /* Input parameters to abso_bg.for */
  double lat, lon;    /* Geographic lat/lon of start-point of ray */
  double freq;        /* radio frequency of ray */
  double elev;        /* forward scattered elevation of ray at ground */
  double R12_index;   /* R12 index */
  int UT[5];          /* UT date and time (YYY, MM, dd, hh, mm) */
  int O_mode, warning_flag;

  /* Output parameters from forward_scatter_loss.for */
  double absorption;

  /* Define the expected input array sizes */
  in_arrsize[0] = 1;
  in_arrsize[1] = 1;
  in_arrsize[2] = 1;
  in_arrsize[3] = 1;
  in_arrsize[4] = 5;
  in_arrsize[5] = 1;
  in_arrsize[6] = 1;

  /* Define the expected output array sizes */
  out_arrsize = 1;

  /* Check for proper number of arguments. */
  if (nrhs != 7) {
    mexErrMsgTxt("ERROR: 7 inputs required when calling abso_bg.");
    return;
  }
  if (nlhs > 1) {
    mexErrMsgTxt("ERROR: 1 output required when calling abso_bg.");
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
  for (ii = 0; ii < 4; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    if (ii == 0) {
      mrows0 = mrows;
      ncols0 = ncols;
    }
    if (mrows != mrows0 || ncols != ncols0) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size does not match first array size.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }
  for (ii = 4; ii < nrhs; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }
  
  /* Check to see if the reference data files exist - abso_bg uses the igrf
     magnetic field routines supplied with IRI2016. */
  if (check_ref_data("iri2016") == 0) return;

  /* Copy the RHS pointers (input from MATLAB) to the input variables of 
     abso_bg_ for those inputs that do not change with lat,lon,elev and freq */
  for (ii = 0; ii < 5; ii++) {
    UT[ii] = (int) *(mxGetPr(prhs[4]) + ii);
  }
  R12_index = mxGetScalar(prhs[5]);
  O_mode = mxGetScalar(prhs[6]);

  /* Create matrix for the return argument. */
  out_arrsize = mrows0 * ncols0;
  plhs[0] = mxCreateDoubleMatrix(mrows0, ncols0, mxREAL);

  for (ray_idx = 0; ray_idx < out_arrsize; ray_idx++) {
    /* Copy the RHS pointers (input from MATLAB) to the input variables of 
       abso_bg_ */
    lat  = *(mxGetPr(prhs[0]) + ray_idx);  
    lon  = *(mxGetPr(prhs[1]) + ray_idx);
    elev = *(mxGetPr(prhs[2]) + ray_idx);
    freq = *(mxGetPr(prhs[3]) + ray_idx);

    /* Call the computational subroutine. */
    absorption = abso_bg_(&lat, &lon, &elev, &freq, UT, &R12_index, &O_mode,
                          &warning_flag);

    if (warning_flag) {
	mexWarnMsgTxt("Absorption potentially unreliable as modified dip > 70 degrees");
    }

    /* Load the output data from abso_bg into the LHS pointers, which are the 
       outputs to MATLAB. */
    *(mxGetPr(plhs[0])+ray_idx) = absorption;
  }
				
}
