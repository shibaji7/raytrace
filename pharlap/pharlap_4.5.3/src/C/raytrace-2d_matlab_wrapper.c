/* 
  =============================================================================
                   raytrace-2d_matlab_wrapper.c
  =============================================================================
   Matlab wrapper for raytrace_2d_
        
   This is a MEX-file for MATLAB.
       
   Change history:
   20/10/05  M.A.Cervera  V1.0  Author
   18/02/10  M.A. Cervera V1.1  Improved error handling
   19/05/11  M.A. Cervera V1.2  Improved efficiency of handling of the 
                                ionospheric and geomagnetic data grids.
   17/01/13  M.A. Cervera V1.3  Now is able to accept the min/max stepsizes
                                in addition to tolerance to control the raytrace
                                ODE solver precision
   17/09/13  M.A. Cervera V1.4  Increased max_pts_in_ray to 20000
   19/02/16  M.A. Cervera V1.5  Extra outputs in ray_path_data
   26/05/16  M.A. Cervera V1.6  Updated for multi-threaded raytrace_2d
   16/05/17  M.A. Cervera V1.7  Now outputs execution time of raytrace call
   12/07/17  M.A. Cervera V1.8  Now doesn't return ray_path_data and ray_state_vec
                                arrays if not requested in the matlab call
   14/09/18  M.A. Cervera V1.9  Now returns the collision frequency and absorption
                                along the ray path and the total absorption.

   13/02/22  M. A. Cervera  V1.10
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings
 
=============================================================================
*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../iono_structures.h"
#include <time.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

void raytrace_2d_(double *origin_lat, double *origin_long, int *num_rays,
		  double *elevs, double *bearing, double *freqs, int *nhops,
		  double *step_size_min, double *step_size_max, double *tol, 
                  struct ionosphere_struct *ionosphere, 
                  int *irregs_flag, int *return_ray_path_data, 
                  int *return_ray_state_vec, double *ray_data, 
                  double *ray_path_data, int *ray_label, int *nhops_attempted,
		  double *ray_state_vec_in, int *npts_in_ray,
		  double *ray_state_vec_out, double *elapsed_time);

void stepmemcpyd(double *dest, double *src, int step, int num_vals);

static struct ionosphere_struct ionosphere;
static int iono_exist_in_mem = 0;

#define max_pts_in_ray 20000    /* maximum number  of points in ray path. If
				   exceeded an error is returned */

/* ----------------------- */
/* The mex gateway routine */
/* ----------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, arrsize, in_arrsize[16], ii, jj, kk, idx, initialize_iono;
  int mrows8, mrows9, mrows10, ncols8, ncols9, ncols10, ncols14, num_bytes;
  int nelem, nfields, new_format, num_numeric;
  int return_ray_path_data, return_ray_state_vec;
  double temp, elapsed_time;
  struct timeval time_start, time_end; 
  char error_mesg[200], filename[200], *refdata_dir;
  char *old_format, *old_format_2D;
  FILE *fopen(), *file_ptr;
  const mwSize *dims, *dims_prev;
  const int ray_data_numfields = 23;
  const char *ray_data_fieldnames[] =
     {"lat", "lon", "ground_range", "group_range", "phase_path",
      "geometric_path_length", "initial_elev", "final_elev", "apogee",
      "gnd_rng_to_apogee", "plasma_freq_at_apogee", "virtual_height",
      "effective_range", "total_absorption", "deviative_absorption", 
      "TEC_path", "Doppler_shift", "Doppler_spread", "FAI_backscatter_loss",
      "frequency", "nhops_attempted", "ray_label", "NRT_elapsed_time"};
  const int ray_data_fieldname_raytrace_output_position[] =
    {0, 1, 2, 3, 16, 10, 6, 7, 4, 5, 13, 15, 11, 18, 12, 17, 9, 8, 14};
  const int ray_path_data_numfields = 11;
  const char *ray_path_data_fieldnames[] =
    {"initial_elev", "frequency", "ground_range", "height", "group_range",
     "phase_path", "geometric_distance", "electron_density", "refractive_index",
     "collision_frequency", "absorption"};
  const int ray_state_vec_numfields = 9;
  const char *ray_state_vec_fieldnames[] =
    {"r", "Q", "theta", "delta_r", "delta_Q", "deviative_absorption",
     "phase_path", "group_path", "group_step_size"};
 
  /* Input parameters to raytrace_2d_ */
  double  origin_lat;           /* geodetic (WGS84) lat. of ray start point */
  double  origin_long;          /* geodetic (WGS84) lon. of ray start point */
  double  *elevs;               /* initial elevation of rays (deg) */
  double  bearing;              /* ray bearing (deg) */
  double  *freqs;               /* carrier frequency of rays (MHz) */
  double  step_size_min;        /* minimum stepsize for ODE solver (km) */ 
  double  step_size_max;        /* maximum stepsize for ODE solver (km) */
  double  tol;                  /* tolerence for each step of raytrace */
  double  *ray_state_vec_in;    /* state vector of the ray starting point */
  int num_rays;                 /* number of input rays */
  int nhops;                    /* number of hops (maximum value: 50) */
  int irregs_flag;              /* flag indicating if doppler spread calculation
                                   is to be performed ( = 0) or not ( = 1) */

  /* Output parameters from raytrace_2d_ */
  double  *ray_data;              /* array containing the ray information */
  double  *ray_path_data;         /* array containing the ray information */
  double  *ray_state_vec_out;     /* state vector of ray at each point 
				     along the ray path */
  int *nhops_attempted;           /* the number of hops actually attempted */
  int *npts_in_ray;               /* the num. points in the ray trajectory */
  int *ray_label;                 /* ray label */


  /* Define the expected input (from matlab) array sizes */
  in_arrsize[0] = 1;       /* latitude */
  in_arrsize[1] = 1;       /* longitude */
  in_arrsize[2] = 1;       /* array of initial elevation for rays (>= 1) */
  in_arrsize[3] = 1;       /* ray bearing */
  in_arrsize[4] = 1;       /* ray carrier frequency for rays (>= 1) */
  in_arrsize[5] = 1;       /* number of hops to perform */
  in_arrsize[6] = 1;       /* 1 or 3, ODE solver tol and min/max stepsize */
  in_arrsize[7] = 1;       /* irregularities flag */
  in_arrsize[8] = max_num_ht*max_num_rng;   /* iono en grid */
  in_arrsize[9] = max_num_ht*max_num_rng;   /* iono en grid 5 min later */
  in_arrsize[10] = max_num_ht*max_num_rng;  /* collision freq. grid */
  in_arrsize[11] = 1;                       /* start height of iono grids */
  in_arrsize[12] = 1;                       /* height increment of iono grids */
  in_arrsize[13] = 1;                       /* range increment of iono grids */
  in_arrsize[14] = 4*max_num_rng;           /* irregularities strength grid */
  in_arrsize[15] = 9;                       /* input ray state vector */

  /* Check for proper number of input arguments. */
  if (nrhs != 8 & nrhs != 9 & nrhs != 15 & nrhs != 16) {
    mexErrMsgTxt("ERROR: incorrect number of input arguments");
    return;
  }

  /* determine if ionosphere need to be initialized */
  initialize_iono = 1;       /*default state - read in the ionosphere */
  if (nrhs == 8 || nrhs == 9) {

    /* OK since we have (nrhs == 8 || nrhs == 9) then an ionosphere has not
       been passed in. Thus, there is not an ionosphere to read in and so we 
       need to set initialize_iono = 0. However, before we do this, check to 
       make sure an ionosphere already exists in memory (i.e. a previous
       raytrace_2d call has already passed it in. If not return an error */
    if (!iono_exist_in_mem) {
      mexErrMsgTxt("ERROR: ionosphere does not exist in memory from a prior raytrace_2d call.");
      return;
    } else {
      initialize_iono = 0;   /* use ionosphere from prior call to raytrace_2d */
    }

  }

  /* Check to ensure required the inputs are numeric. The optional input is
     checked later. */
  if (nrhs == 8 || nrhs == 15) {
    num_numeric = nrhs;
  } else {
    num_numeric = nrhs-1;    /* for this case the last input is the optional */
  }                          /* state vecor of the rays which is a structure */
  for (ii = 0; ii < num_numeric; ii++) {
    if (mxIsNumeric(prhs[ii]) == 0) {
      snprintf(error_mesg, 200, "ERROR: Input %d must be numeric.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }
  
  /* Check to ensure the optional input (if supplied) is a structure and has */
  /* the correct size  */
  if (nrhs == 9 || nrhs == 16) {
    ii = nrhs - 1;
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    dims = mxGetDimensions(prhs[2]);   
    num_rays = dims[1];     /* this is the number of input rays */
    if(!mxIsStruct(prhs[ii])) {
       mexErrMsgIdAndTxt( "MATLAB:ray_state_vec_in:inputNotStruct",
			  "ERROR: Input %d must be a structure.", ii+1);
       return;
    }
    nfields = mxGetNumberOfFields(prhs[ii]);
    if (nfields != 9) {
      snprintf(error_mesg, 200, 
	  "ERROR: Input %d has incorrect number of fields in structure.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (ncols != num_rays) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* Get and check the size of the input scalers and arrays. */

  /* make sure inputs 0 to 7 have the correct size */
  for (ii = 0; ii < 8; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (ii == 2 & arrsize != 0) arrsize = 1;  /* must be >= 1 */
    if (ii == 4 & arrsize != 0) arrsize = 1;  /* must be >= 1 */
    if (ii == 6 & arrsize == 3) arrsize = 1;  /* may be 1 or 3 elements */
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }  

  /* make sure inputs 2 and 4 have the same size */
  dims_prev = mxGetDimensions(prhs[2]);
  dims = mxGetDimensions(prhs[4]);
  if (dims[0] != dims_prev[0] || dims[1] != dims_prev[1]) {
    snprintf(error_mesg, 200,
	    "ERROR: Array size of inputs 3, and 5 are inconsistent.");
    mexErrMsgTxt(error_mesg);
    return;
  }

  /* if the ionospheric grids have been passed in then check their validity */
  if (initialize_iono) {

    /* make sure inputs 8, 9 and 10 are not too large */
    for (ii = 8; ii < 11; ii++) { 
      mrows = mxGetM(prhs[ii]);
      ncols = mxGetN(prhs[ii]);
      if (mrows > max_num_ht || ncols > max_num_rng) {
	snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
	mexErrMsgTxt(error_mesg);
	return;
      }
    }

    /* make sure inputs 11, 12 and 13 have the correct size */
    for (ii = 11; ii < 14; ii++) {
      mrows = mxGetM(prhs[ii]);
      ncols = mxGetN(prhs[ii]);
      arrsize = mrows*ncols;
      if (arrsize != in_arrsize[ii]) {
	snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
	mexErrMsgTxt(error_mesg);
	return;
      }
    }

    /* make sure input 14 has the correct size */
    ii = 14;
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    if (mrows != 4 || ncols > max_num_rng) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }

    /* make sure array size of inputs 8, 9, 10 and 14 are consistent */
    mrows8 = mxGetN(prhs[8]);
    mrows9 = mxGetN(prhs[9]);
    mrows10 = mxGetN(prhs[10]);
    if (mrows9 - mrows8 != 0 || mrows10 - mrows8 != 0) {
      snprintf(error_mesg, 200, "ERROR: Number of rows of inputs 9, 10 and 11 are inconsistent.");
      mexErrMsgTxt(error_mesg);
      return;
    }

    ncols8 = mxGetN(prhs[8]);
    ncols9 = mxGetN(prhs[9]);
    ncols10 = mxGetN(prhs[10]);
    ncols14 = mxGetN(prhs[14]);
    if (ncols9 - ncols8 != 0 || ncols10 - ncols8 != 0 || ncols14 - ncols8 != 0) {
      snprintf(error_mesg, 200, "ERROR: Number of columns of inputs 9, 10, 11, and 15 are inconsistent.");
      mexErrMsgTxt(error_mesg);
      return;
    }

  }

  /* get the number of hops and make sure that it is not too large and > 0 */
  nhops = mxGetScalar(prhs[5]);
  if (nhops > 50) {
    mexErrMsgTxt("ERROR: number of hops is too large, it must be < 50");
    return;
  }
  if (nhops < 1) {
    mexErrMsgTxt("ERROR: number of hops must be >= 1");
    return;
  }
  
  /* Copy the RHS pointers (input from MATLAB) to the input variables of 
     raytrace_2d_ */
  origin_lat = mxGetScalar(prhs[0]);
  origin_long = mxGetScalar(prhs[1]);
  dims = mxGetDimensions(prhs[2]);
  num_rays = dims[1];
  elevs = mxMalloc(num_rays * sizeof(double));
  for (ii = 0; ii < num_rays; ii++) {
     elevs[ii] = *(mxGetPr(prhs[2]) + ii);
  }  
  bearing = mxGetScalar(prhs[3]);
  freqs = mxMalloc(num_rays * sizeof(double));
  for (ii = 0; ii < num_rays; ii++) {
    freqs[ii] = *(mxGetPr(prhs[4]) + ii);
  }
  nhops = mxGetScalar(prhs[5]);
  irregs_flag = mxGetScalar(prhs[7]);

  mrows = mxGetM(prhs[6]);
  ncols = mxGetN(prhs[6]);
  nelem = mrows*ncols;
  if (nelem == 1) {
    temp = mxGetScalar(prhs[6]);
    if (temp > 1e-12 & temp < 1e-2) {
      tol = temp;
      step_size_min = (double) 0.01;
      step_size_max = (double) 10.0;
    }
    else if (temp == 1) {
      tol = 1e-8;
      step_size_min = (double) 0.01;
      step_size_max = (double) 10.0;
    }
    else if (temp == 2) {
      tol = 1e-7;
      step_size_min = (double) 0.025;
      step_size_max = (double) 25.0;
    }
    else if (temp == 3) {
      tol = 1e-6;
      step_size_min = (double) 0.1;
      step_size_max = (double) 100.0;
    }
    else {
      mexErrMsgTxt("ERROR: Incorrect values for min/max stepsize and tolerance");
      return;
    }      
  }
  else {
    tol = *(mxGetPr(prhs[6]) + 0);
    step_size_min = *(mxGetPr(prhs[6]) + 1);
    step_size_max = *(mxGetPr(prhs[6]) + 2);
  }
  if (tol < 1e-12 || tol > 1e-2 || step_size_min > step_size_max ||
      step_size_min < 0.001 || step_size_min > 1 || step_size_max > 100 ||
      step_size_max < 1) {
    mexErrMsgTxt("ERROR: Incorrect values for min/max stepsize and tolerance");
    return;
  }           

  /* if the ionosphere is required to be initialized (rather than using that 
     from a prior call to raytrace_2d) then read it in from Matlab */
  if (initialize_iono) {
    mrows = mxGetM(prhs[8]);
    ncols = mxGetN(prhs[8]);
    for (jj = 0; jj < ncols; jj++) {
      for (ii = 0; ii < mrows; ii++) {
	idx = jj*mrows + ii;
	ionosphere.eN[jj][ii] = *(mxGetPr(prhs[8]) + idx);
	ionosphere.eN_5[jj][ii] = *(mxGetPr(prhs[9]) + idx);
	ionosphere.col_freq[jj][ii] = *(mxGetPr(prhs[10]) + idx);
      }
    }

    mrows = mxGetM(prhs[14]);
    ncols = mxGetN(prhs[14]);
    for (jj = 0; jj < ncols; jj++) {
      ionosphere.irreg_strength[jj] = *(mxGetPr(prhs[14]) + jj*mrows);
      ionosphere.irreg_sma_dip[jj] = *(mxGetPr(prhs[14]) + jj*mrows + 1);
      ionosphere.irreg_sma_azim[jj] = *(mxGetPr(prhs[14]) + jj*mrows + 2);
      ionosphere.dop_spread_sq[jj] = *(mxGetPr(prhs[14]) + jj*mrows+3);
    } 

    ionosphere.nRange = mxGetN(prhs[9]);
    ionosphere.NumHt = mxGetM(prhs[9]);
    ionosphere.HtMin = mxGetScalar(prhs[11]);
    ionosphere.HtInc = mxGetScalar(prhs[12]);
    ionosphere.dRange = mxGetScalar(prhs[13]);

    /* Now the ionosphere has been read in set the iono_exist_in_mem flag to
       indicate this for future raytrace calls */
    iono_exist_in_mem = 1;

  }


  /* read in the optional input (structure containing a user defined starting 
     ray state vector for each ray) from  Matlab (if required) and check to make
     sure each field is valid */
  mxArray *tmp_ptr;
  ray_state_vec_in = (double *) mxMalloc(9*num_rays*sizeof(double));
  if (nrhs == 9 || nrhs == 16) {
    for (kk = 0; kk < num_rays; kk++) {
      for (jj = 0; jj < 9; jj++) {
	idx = kk*9 + jj;
	tmp_ptr = mxGetFieldByNumber(prhs[nrhs-1],(mwIndex) kk, (mwIndex) jj);
	if (tmp_ptr == NULL) {
	  mexPrintf("\n%s%d\t%s%d\t%s%d\n", "Input: ", nrhs, 
		    "Structure Index: ", kk+1, "Field: ", jj+1);
	  mexErrMsgIdAndTxt( "MATLAB:ray_state_vec_in:fieldEmpty",
			     "Above field is empty!");
	  return;
	}
	if (!mxIsNumeric(tmp_ptr) || mxIsComplex(tmp_ptr)) {
	  mexPrintf("\n%s%d\t%s%d\t%s%d\n", "Input: ", nrhs, 
		    "Structure Index :", kk+1, "Field: ", jj+1);
	  mexErrMsgIdAndTxt( "MATLAB:ray_state_vec_in:fieldInvalid",
			     "Above field is not a valid numeric value!");
	  return;
	}
	memcpy(ray_state_vec_in + idx, mxGetData(tmp_ptr), sizeof(double));
	if (mxIsNaN(*(ray_state_vec_in +idx)) |
	    mxIsInf(*(ray_state_vec_in +idx))) {
	  mexPrintf("\n%s%d\t%s%d\t%s%d\n", "Input: ", nrhs, 
		    "Structure Index :", kk+1, "Field: ", jj+1);
	  mexErrMsgIdAndTxt( "MATLAB:ray_state_vec_in:fieldInvalid",
			     "Above field is not a valid numeric value!");
	return;

	}
      }
    } 
  }
  else {   /* ray_state_vec_in not input from Matlab so set all values to -1 */
    for (kk = 0; kk < num_rays; kk++) {
      for (jj = 0; jj < 9; jj++) {
	idx = kk*9 + jj;
        *(ray_state_vec_in + idx) = -1;
      }
    }
  }

  /* malloc space for the raytracing output */
  nhops_attempted = (int *) mxMalloc(num_rays*sizeof(int));
  npts_in_ray = (int *) mxMalloc(num_rays*sizeof(int));
  ray_data = (double *) mxMalloc(19*nhops*num_rays*sizeof(double));
  ray_path_data =(double *) mxMalloc(9*max_pts_in_ray*num_rays*sizeof(double));
  ray_label = (int *) mxMalloc(nhops*num_rays*sizeof(int));
  ray_state_vec_out =
                (double *) mxMalloc(9*max_pts_in_ray*num_rays*sizeof(double));
  
  /* determine if the ray_path_data and ray_state_vec arrays have been 
     requested to be returned to matlab */
  return_ray_path_data = 0;      /* default action - don't return this array */
  return_ray_state_vec = 0;      /* default action - don't return this array */
  if (nlhs == 2) {
    return_ray_path_data = 1;    /* return the ray_path_data array */
  }
  if (nlhs == 3) {
    return_ray_path_data = 1;    /* return the ray_path_data array */
    return_ray_state_vec = 1;    /* return the ray_state_vec array */
  }

  /* Call the computational subroutine. */
  raytrace_2d_(&origin_lat, &origin_long, &num_rays, elevs, &bearing, freqs,
	       &nhops, &step_size_min, &step_size_max, &tol, &ionosphere,
	       &irregs_flag, &return_ray_path_data, &return_ray_state_vec, 
               ray_data, ray_path_data, ray_label, nhops_attempted, 
               ray_state_vec_in, npts_in_ray, ray_state_vec_out, &elapsed_time);

  /* Copy the raytrace data into the matlab data structures. */
  /* If only 1 ray has been raytraced (num_rays == 1) and either of the 
     environment variables PHARLAP_OLD_FORMAT or PHARLAP_OLD_FORMAT_2D has been 
     set to "true" then return the data in Matlab Double Matricies. This is to 
     retain backwards compatibility with Matlab code using PHaRLAP version 
     3.7.1 and earlier. */
  new_format = 1;               /* new format is the default */
  old_format = getenv("PHARLAP_OLD_FORMAT");
  old_format_2D = getenv("PHARLAP_OLD_FORMAT_2D");
  
  if (old_format != NULL) {
    if (strcmp(old_format, "true") == 0 && num_rays == 1) new_format = 0;
  }
  if (old_format_2D != NULL) {
    if (strcmp(old_format_2D, "true") == 0 && num_rays == 1) new_format = 0;
  }
  
  if (new_format) {     /* copy data to Matlab structures */

    /* Create the various structures for the return arguments to matlab */ 
    /* ray data */
    plhs[0] = mxCreateStructMatrix((mwSize) 1, (mwSize) num_rays,
				   ray_data_numfields, ray_data_fieldnames);

    /* ray path data*/
    plhs[1] = mxCreateStructMatrix((mwSize) 1, (mwSize) num_rays,
			    ray_path_data_numfields, ray_path_data_fieldnames);

    /* ray state vector data */
    plhs[2] = mxCreateStructMatrix((mwSize) 1, (mwSize) num_rays,
			   ray_state_vec_numfields, ray_state_vec_fieldnames); 

    /* loop over rays and copy data to Matlab structures */
    tmp_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], (mwIndex) 0, "NRT_elapsed_time", tmp_ptr);
    *mxGetPr(tmp_ptr) = elapsed_time;
    for (kk = 0; kk < num_rays; kk++) {

      /* copy the ray_data structure */
      for (jj = 0; jj < ray_data_numfields-4; jj++) {
	idx = kk + num_rays* ray_data_fieldname_raytrace_output_position[jj];
	tmp_ptr = mxCreateDoubleMatrix(1, nhops_attempted[kk], mxREAL);
	mxSetFieldByNumber(plhs[0], (mwIndex) kk, (mwIndex) jj, tmp_ptr);
	stepmemcpyd(mxGetPr(tmp_ptr), &ray_data[idx], num_rays*19, 
	                                                  nhops_attempted[kk]);
      }
      tmp_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
      mxSetField(plhs[0], (mwIndex) kk, "frequency", tmp_ptr);
      *mxGetPr(tmp_ptr) = freqs[kk];
      
      tmp_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
      mxSetField(plhs[0], (mwIndex) kk, "nhops_attempted", tmp_ptr);
      *mxGetPr(tmp_ptr) = nhops_attempted[kk];

      tmp_ptr = mxCreateDoubleMatrix(1, nhops_attempted[kk], mxREAL);
      mxSetField(plhs[0], (mwIndex) kk, "ray_label", tmp_ptr);
      for (ii = 0; ii < nhops_attempted[kk]; ii++) {
	idx = kk + num_rays*ii;
	*(mxGetPr(tmp_ptr) + ii) = (double) ray_label[idx];
      }

      /* copy the ray_path_data structure if requested */
      if (return_ray_path_data) {
        tmp_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[1], (mwIndex) kk, "initial_elev", tmp_ptr);
	*mxGetPr(tmp_ptr) = elevs[kk];

	tmp_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[1], (mwIndex) kk, "frequency", tmp_ptr);
	*mxGetPr(tmp_ptr) = freqs[kk];

	for (jj = 0; jj < ray_path_data_numfields-2; jj++) {
	  idx = kk + num_rays*jj;
	  tmp_ptr = mxCreateDoubleMatrix(1, npts_in_ray[kk], mxREAL);
	  mxSetFieldByNumber(plhs[1], (mwIndex) kk, (mwIndex) jj+2, tmp_ptr);
	  stepmemcpyd(mxGetPr(tmp_ptr), &ray_path_data[idx], num_rays*9, 
							       npts_in_ray[kk]);
	}
      }

      /* copy the ray_state_vec structure if requested */
      if (return_ray_state_vec) {
        for (jj = 0; jj < ray_state_vec_numfields; jj++) {
	  idx = kk + num_rays*jj;
	  tmp_ptr = mxCreateDoubleMatrix(1, npts_in_ray[kk], mxREAL);
	  mxSetFieldByNumber(plhs[2], (mwIndex) kk, (mwIndex) jj, tmp_ptr);
	  stepmemcpyd(mxGetPr(tmp_ptr), &ray_state_vec_out[idx], num_rays*9,
		                                           npts_in_ray[kk]);
        }
      }
      
    }
    
  } else {          /* copy data to Matlab Double Matricies */

    /* Create matricies for the return arguments. */
    plhs[0] = mxCreateDoubleMatrix(19, nhops_attempted[0], mxREAL);
    plhs[1] = mxCreateDoubleMatrix(9, npts_in_ray[0], mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, nhops_attempted[0], mxREAL);
    plhs[4] = mxCreateDoubleMatrix(9, npts_in_ray[0], mxREAL);

    /* Load the output data from raytrace into the LHS pointers, which are the 
       outputs to MATLAB. */
    for (ii=0; ii< nhops_attempted[0]; ii++) {
      idx = ii*19;
      stepmemcpyd(mxGetPr(plhs[0])+idx, &ray_data[ii], nhops, 19);
    }
    for (ii=0; ii< npts_in_ray[0]; ii++) {
      idx = ii*9;
      stepmemcpyd(mxGetPr(plhs[1])+idx, &ray_path_data[ii], max_pts_in_ray, 9);
    }
    *mxGetPr(plhs[2]) = nhops_attempted[0];

    for (jj = 0; jj < nhops_attempted[0]; jj++) {
      *(mxGetPr(plhs[3]) + jj) = (double) (&ray_label[0])[jj];
    }

    for (ii=0; ii< npts_in_ray[0]; ii++) {
      idx = ii*9;
      stepmemcpyd(mxGetPr(plhs[4])+idx, &ray_state_vec_out[ii],
						      max_pts_in_ray, 9);
    }   
    
  }

  /* cleanup - free malloced space */
  mxFree(elevs);
  mxFree(freqs);
  mxFree(nhops_attempted);
  mxFree(npts_in_ray);
  mxFree(ray_state_vec_in);
  mxFree(ray_data);
  mxFree(ray_path_data);
  mxFree(ray_state_vec_out);
  mxFree(ray_label);
  
}

