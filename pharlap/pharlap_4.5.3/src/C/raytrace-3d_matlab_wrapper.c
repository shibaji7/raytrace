/* 
  =============================================================================
                   raytrace-3d_matlab_wrapper.c
  ==============================================================================
   Matlab wrapper for raytrace_3d_
        
   This is a MEX-file for MATLAB.
       
   Change history:
   07/12/09  M.A. Cervera V1.0  Author
   18/02/10  M.A. Cervera V1.1  Improved error handling
   19/05/11  M.A. Cervera V1.2  Improved efficiency of handling of the 
                                ionospheric and geomagnetic data grids.
   17/01/13  M.A. Cervera V1.3  Now is able to accept the min/max stepsizes
                                in addition to tolerance to control the raytrace
                                ODE solver precision
   04/08/13  M.A. Cervera V1.4  Increased max_pts_in_ray to 20000 and fixed
                                bug with checking for valid input arguments.
   24/07/14  M.A. Cervera V1.5  Increased the ionospheric grid sizes from 
                                201x201x401 to 701x701x301 lat, lon, heights
   22/04/15  M.A. Cervera V1.6  Increased the ionospheric grid sizes from 
                                701x701x401 to 701x701x401 lat, lon, heights
   28/04/15  M.A. Cervera V1.7  Added routine called from mexAtExit to clear
                                ionosphere when clear is called from Matlab
   03/12/15  M.A. Cervera V2.0  Now calls multi-threaded raytrace_3d
   23/05/16  M.A. Cervera V2.1  Fixed minor bug where the last required input
                                was not being properly checked 
   10/06/16  M.A. Cervera V2.2  Additional fields added to outputs
   14/10/16  M.A. Cervera V2.3  Additional input checking
   16/05/17  M.A. Cervera V2.4  Now outputs execution time of raytrace call
   12/07/17  M.A. Cervera V2.5  Now doesn't return ray_path_data and 
                                ray_state_vec  arrays if not requested in the 
                                matlab call
   25/01/18  M.A. Cervera V2.6  Fixed a memory leak issue when freeing allocated
                                memeory for the ionosphere via the 
                                clear_ionosphere() function
   05/06/18  M.A. Cervera V2.7  Now returns the collision frequency along the
                                ray path
   10/07/18  M.A. Cervera V2.8  Now returns the total absorption 

   13/02/22  M. A. Cervera  V2.9
     Now using snprintf (under unix) and sprintf_s (under Windows)
     instead of sprintf for generating strings

  ==============================================================================
*/

#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../iono_structures_3d.h"
#include <time.h>

#if defined(WIN64)
  #define snprintf sprintf_s    /* needed for gcc under Windows/Cygwin */
#endif

void raytrace_3d_(double *start_lat, double *start_long, double *start_height,
                  int *num_rays, double *elevs, double *bearings,
		  double *freqs, int *OX_mode, int *nhops,
		  double *step_size_min, double *step_size_max, double *tol,
		  struct ionosphere_struct *ionosphere,
                  struct geomag_field_struct *geomag_field,  
                  double *ray_state_vec_in, int *return_ray_path_data, 
	          int *return_ray_state_vec, double *ray_data, 
                  double *ray_path_data, int *ray_label, int *nhops_attempted, 
                  int *npts_in_ray, double *ray_state_vec_out, 
                  double *elapsed_time);

void stepmemcpyd(double *dest, double *src, int step, int num_vals);

/* static struct ionosphere_struct ionosphere;  - ionospheric data too large and
   so this statement was causing the stack to be "squashed" ie corrupted - so we
   want to allocate from heap instead */
static struct ionosphere_struct *ptr_ionosphere = NULL;
static struct geomag_field_struct geomag_field;
static int iono_exist_in_mem = 0;

/* clear the ionosphere from memory when mexAtExit is triggered by a clear 
   command at the Matlab prompt */
static void clear_ionosphere(void)
{
  if (ptr_ionosphere != NULL)
  {
    mxFree(ptr_ionosphere);
    ptr_ionosphere = NULL;
  }
}

#define max_pts_in_ray 20000    /* maximum number  of points in ray path. If
				   exceeded an error is returned */

/* ----------------------- */
/* The mex gateway routine */
/* ----------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mrows, ncols, arrsize, in_arrsize[23], ii, jj, kk, idx, nelem;
  int mrows8, mrows9, mrows10, ncols8, ncols9, ncols10, ncols14;
  int initialize_iono, iono_grid_size, nfields, new_format, num_numeric;
  int return_ray_path_data, return_ray_state_vec;
  double temp, elapsed_time;
  struct timeval time_start, time_end; 
  char error_mesg[200], filename[200], *refdata_dir;
  char *old_format, *old_format_3D;
  FILE *fopen(), *file_ptr;
  const mwSize *dims, *dims_prev;
  mwSize number_of_dimensions;
  const int ray_data_numfields = 19;
  const char *ray_data_fieldnames[] =
     {"lat", "lon", "ground_range", "group_range", "phase_path",
      "initial_elev", "final_elev", "initial_bearing", "final_bearing",
      "total_absorption", "deviative_absorption", "TEC_path", "Doppler_shift", 
      "apogee", "geometric_path_length", "frequency", "nhops_attempted", 
      "ray_label", "NRT_elapsed_time"};
  const int ray_data_fieldname_raytrace_output_position[] =
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  const int ray_path_data_numfields = 22;
  const char *ray_path_data_fieldnames[] =
    {"initial_elev", "initial_bearing", "frequency", "lat", "lon", "height",
     "group_range", "phase_path", "refractive_index", "group_refractive_index",
     "wavenorm_ray_angle", "wavenorm_B_angle", "polariz_mag",
     "wave_Efield_tilt", "volume_polariz_tilt", "electron_density", "geomag_x",
     "geomag_y", "geomag_z", "geometric_distance", "collision_frequency", 
     "absorption"};  
  const int ray_state_vec_numfields = 12;
  const char *ray_state_vec_fieldnames[] =
    {"pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "group_path",
     "geometrical_path", "phase_path", "absorption", "indep_var", "ODE_step_size"};
  
  /* Input parameters to raytrace_3d_ */
  double  start_lat;            /* geodetic (WGS84) lat. of ray start point */
  double  start_long;           /* geodetic (WGS84) lon. of ray start point */
  double  start_height;         /* geodetic (WGS84) height of ray start point */
  double  *elevs;               /* initial elevation of rays (deg) */
  double  *bearings;            /* initial ray bearings (deg) */
  double  *freqs;               /* ray carrier frequency (MHz) */
  double  step_size_min;        /* minimum stepsize for ODE solver (km) */ 
  double  step_size_max;        /* maximum stepsize for ODE solver (km) */
  double  tol;                  /* tolerence for each step of raytrace */
  double  *ray_state_vec_in;    /* state vector at the rays' starting point */
  int nhops;                    /* number of hops (maximum value: 50) */
  int OX_mode;                  /* polarization mode of ray: 1 = O, -1 = X, 
                                   0 = no field */
  int num_rays;

  /* Output parameters from raytrace_3d_ */
  double  *ray_data;              /* array containing the ray information */
  double  *ray_path_data;         /* array containing the ray information */
  double  *ray_state_vec_out;     /* state vector of ray at each point 
				     along the ray path */
  int *nhops_attempted;           /* the number of hops actually attempted */
  int *npts_in_ray;               /* the num. points in the ray trajectory */
  int *ray_label;                 /* ray label */
 
  iono_grid_size = max_num_ht * max_num_lat * max_num_lon;
  
  /* Define the expected input (from matlab) array sizes */
  in_arrsize[0] = 1;               /* start latitude */
  in_arrsize[1] = 1;               /* start longitude */
  in_arrsize[2] = 1;               /* start height */
  in_arrsize[3] = 1;               /* minimum number of start elevations */
  in_arrsize[4] = 1;               /* minimum number of ray bearings */
  in_arrsize[5] = 1;               /* minimum number of frequencies */
  in_arrsize[6] = 1;               /* OX_mode */
  in_arrsize[7] = 1;               /* number of hops to complete */
  in_arrsize[8] = 1;               /* ODE solver tolerance - can also be 3 */
  in_arrsize[9] = iono_grid_size;  /* ionospheric elec density grid */
  in_arrsize[10] = iono_grid_size; /* elec density grid 5 mins later */
  in_arrsize[11] = iono_grid_size; /* ionospheric collision frequency grid */
  in_arrsize[12] = 9;              /* ionospheric grid parameters */
  in_arrsize[13] = 101*101*201;    /* gridded Bx */
  in_arrsize[14] = 101*101*201;    /* gridded By */
  in_arrsize[15] = 101*101*201;    /* gridded Bz */
  in_arrsize[16] = 9;              /* geomagnetic grid parameters */
  in_arrsize[17] = 12;             /* optional input - state vector at start */
  
  /* Check for proper number of input arguments. */
  if (nrhs != 9 && nrhs != 10 && nrhs != 17 && nrhs != 18) {
    mexErrMsgTxt("ERROR: incorrect number of input arguments");
    return;
  }

  /* determine if ionospheric and geomagnetic grids need to be initialized */
  if (ptr_ionosphere == NULL) {
    mexAtExit(clear_ionosphere);
    ptr_ionosphere = mxMalloc(sizeof(*ptr_ionosphere));
    mexMakeMemoryPersistent(ptr_ionosphere);
  }
  initialize_iono = 1;     /*default state - read in the ionosphere */
  if (nrhs == 9 || nrhs == 10) {

    /* OK since we have (nrhs == 9 || nrhs == 10) then an ionosphere has not
       been passed in. Thus, there is not an ionosphere to read in and so we 
       need to set initialize_iono = 0. However, before we do this, check to 
       make sure an ionosphere already exists in memory (i.e. a previous
       raytrace_3d call has already passed it in. If not return an error */
    if (!iono_exist_in_mem) {
      mexErrMsgTxt("ERROR: ionosphere does not exist in memory from a prior raytrace_3d call.");
      return;
    } else {
      initialize_iono = 0;   /* use ionospheric and geomagnetic data grids from 
                                prior call to raytrace_3d */
    }
  }

  /* Check to ensure the required inputs are numeric. The optional input is 
     checked later.*/
  if (nrhs == 9 || nrhs == 17) {
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
  if (nrhs == 10 || nrhs == 18) {
    ii = nrhs - 1;
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    dims = mxGetDimensions(prhs[3]);   
    num_rays = dims[1];     /* this is the number of input rays */
    if(!mxIsStruct(prhs[ii])) {
       mexErrMsgIdAndTxt( "MATLAB:ray_state_vec_in:inputNotStruct",
			  "ERROR: Input %d must be a structure.", ii+1);
       return;
    }
    nfields = mxGetNumberOfFields(prhs[ii]);
    if (nfields != 12) {
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

  /* make sure inputs 0 to 2 have the correct size */
  for (ii = 0; ii < 3; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }
  
  /* make sure array size of inputs 3 to 5 are consistent and not empty*/
  if (mxIsEmpty(prhs[3])) { 
    snprintf(error_mesg, 200, "ERROR: input 4 is an empty array");
    mexErrMsgTxt(error_mesg);
    return;
  }
  dims_prev = mxGetDimensions(prhs[3]);
  for (ii = 4; ii < 6; ii++) {
    dims = mxGetDimensions(prhs[ii]);
    if (dims[0] != dims_prev[0] || dims[1] != dims_prev[1]) {
      snprintf(error_mesg, 200,
	      "ERROR: Array size of inputs 4, 5 and 6 are inconsistent.");
      mexErrMsgTxt(error_mesg);
      return;
    }
    dims_prev = mxGetDimensions(prhs[ii]);
  }
  
  /* make sure inputs 6 to 8 have the correct size */
  for (ii = 6; ii < 9; ii++) {
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (ii == 8 && arrsize == 3) arrsize = 1;  /* may be 1 or 3 elements */
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* if the ionospheric and geomagnetic grids have been passed in then check 
     their validity */
  if (initialize_iono) {

    /* make sure inputs 9 to 11 have the correct number of dimensions and
       are not too large */ 
    for (ii = 9; ii < 12; ii++) { 
      dims = mxGetDimensions(prhs[ii]);
      number_of_dimensions = mxGetNumberOfDimensions(prhs[ii]);
      if (number_of_dimensions != 3 || dims[0] > max_num_lat ||
	  dims[1] > max_num_lon || dims[2] > max_num_ht) {
	snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
	mexErrMsgTxt(error_mesg);
	return;
      }
    }

    /* make sure array size of inputs 9 to 11 are consistent */ 
    dims_prev = mxGetDimensions(prhs[9]);
    for (ii = 10; ii < 12; ii++) { 
      dims = mxGetDimensions(prhs[ii]);
      if (dims[0] != dims_prev[0] || dims[1] != dims_prev[1] || 
	  dims[2] != dims_prev[2]) {
	snprintf(error_mesg, 200,
		"ERROR: Array size of inputs 10, 11 and 12 are inconsistent.");
	mexErrMsgTxt(error_mesg);
	return;
      }
      dims_prev = mxGetDimensions(prhs[ii]);
    }

    /* make sure input 12 has the correct size */
    ii = 12;
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }

    /* make sure inputs 13 to 15 have the correct number of dimensions and
       are not too large */ 
    for (ii = 13; ii < 15; ii++) { 
      dims = mxGetDimensions(prhs[ii]);
      number_of_dimensions = mxGetNumberOfDimensions(prhs[ii]);
      if (number_of_dimensions != 3 || dims[0] > 101 || dims[1] > 101 || 
	  dims[2] > 201) {
	snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
	mexErrMsgTxt(error_mesg);
	return;
      }
    }

    /* make sure array size of inputs 13 to 15 are consistent */ 
    dims_prev = mxGetDimensions(prhs[13]);
    for (ii = 14; ii < 16; ii++) { 
      dims = mxGetDimensions(prhs[ii]);
      if (dims[0] != dims_prev[0] || dims[1] != dims_prev[1] || 
	  dims[2] != dims_prev[2]) {
	snprintf(error_mesg, 200,
		"ERROR: Array size of inputs 14, 15 and 16 are inconsistent.");
	mexErrMsgTxt(error_mesg);
	return;
      }
      dims_prev = mxGetDimensions(prhs[ii]);
    }

    /* make sure input 16 has the correct size */
    ii = 16;
    mrows = mxGetM(prhs[ii]);
    ncols = mxGetN(prhs[ii]);
    arrsize = mrows*ncols;
    if (arrsize != in_arrsize[ii]) {
      snprintf(error_mesg, 200, "ERROR: Input %d array size is incorrect.", ii+1);
      mexErrMsgTxt(error_mesg);
      return;
    }

  }
  
  /* get the number of hops and make sure that it is not too large and > 0 */
  nhops = mxGetScalar(prhs[7]);
  if (nhops > 50) {
    mexErrMsgTxt("ERROR: number of hops is too large, it must be < 50");
    return;
  }
  if (nhops < 1) {
    mexErrMsgTxt("ERROR: number of hops must be >= 1");
    return;
  }
    
  /* Copy the RHS pointers (input from MATLAB) to the input variables of 
     raytrace_3d_ */
  start_lat = mxGetScalar(prhs[0]);
  if (start_lat < -90.0 || start_lat > 90.0) {
      mexErrMsgTxt("The start latitude of the ray must be in the range -90 to 90 degrees.");
      return;
    }
  start_long = mxGetScalar(prhs[1]);
  if (start_long < -180.0 || start_long > 180.0) {
      mexErrMsgTxt("The start longitude of the ray must be in the range -180 to 180 degrees.");
      return;
    }
  start_height = mxGetScalar(prhs[2]);
  dims = mxGetDimensions(prhs[3]);
  num_rays = dims[1];
  elevs = mxMalloc(num_rays * sizeof(double));
  bearings = mxMalloc(num_rays * sizeof(double));
  freqs = mxMalloc(num_rays * sizeof(double));
  for (ii = 0; ii < num_rays; ii++) {
     elevs[ii] = *(mxGetPr(prhs[3]) + ii);
     bearings[ii] = *(mxGetPr(prhs[4]) + ii);
     freqs[ii] = *(mxGetPr(prhs[5]) + ii);
  }
  OX_mode = mxGetScalar(prhs[6]);
  nhops = mxGetScalar(prhs[7]);
  mrows = mxGetM(prhs[8]);
  ncols = mxGetN(prhs[8]);
  nelem = mrows*ncols;
  if (nelem == 1) {
    temp = mxGetScalar(prhs[8]);
    if (temp > 1e-12 && temp < 1e-2) {
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
    tol = *(mxGetPr(prhs[8]) + 0);
    step_size_min = *(mxGetPr(prhs[8]) + 1);
    step_size_max = *(mxGetPr(prhs[8]) + 2);
  }
  if (tol < 1e-12 || tol > 1e-2 || step_size_min > step_size_max ||
      step_size_min < 0.001 || step_size_min > 1 || step_size_max > 100 ||
      step_size_max < 1) {
    mexErrMsgTxt("ERROR: Incorrect values for min/max stepsize and tolerance");
    return;
  }             

  /* if the ionosphere is required to be initialized (rather than using that 
     from a prior call to raytrace_3d) then read it in from Matlab */
  if (initialize_iono) {
    dims = mxGetDimensions(prhs[9]);
    for (kk = 0; kk < dims[2]; kk++) {
      for (jj = 0; jj < dims[1]; jj++) {
	for (ii = 0; ii < dims[0]; ii++) {
	  idx = kk*dims[1]*dims[0] + jj*dims[0] + ii;
	  ptr_ionosphere->eN[kk][jj][ii] = *(mxGetPr(prhs[9]) + idx);
	  ptr_ionosphere->eN_5[kk][jj][ii] = *(mxGetPr(prhs[10]) + idx);
	  ptr_ionosphere->col_freq[kk][jj][ii] = *(mxGetPr(prhs[11]) + idx);
	}
      }
    }
    ptr_ionosphere->lat_min = *(mxGetPr(prhs[12]) + 0);
    if (ptr_ionosphere->lat_min < -90.0 || ptr_ionosphere->lat_min > 90.0) {
      mexErrMsgTxt("The start latitude of the input ionospheric grid must be in the range -90 to 90 degrees.");
      return;
    }
    ptr_ionosphere->lat_inc = *(mxGetPr(prhs[12]) + 1);
    ptr_ionosphere->num_lat = *(mxGetPr(prhs[12]) + 2);
    if (dims[0] != ptr_ionosphere->num_lat) {
      mexErrMsgTxt("The number of latitudes in the input ionospheric grid does not match the input grid parameter.");
      return;
    }
    ptr_ionosphere->lat_max = ptr_ionosphere->lat_min + 
		       (ptr_ionosphere->num_lat - 1) * ptr_ionosphere->lat_inc;
    ptr_ionosphere->lon_min = *(mxGetPr(prhs[12]) + 3);
    if (ptr_ionosphere->lon_min < -180.0 || ptr_ionosphere->lon_min > 180.0) {
      mexErrMsgTxt("The start longitude of the input ionospheric grid must be in the range -180 to 180 degrees.");
      return;
    }
    ptr_ionosphere->lon_inc = *(mxGetPr(prhs[12]) + 4);
    ptr_ionosphere->num_lon = *(mxGetPr(prhs[12]) + 5);
    if (dims[1] != ptr_ionosphere->num_lon) {
      mexErrMsgTxt("The number of longitudes in the input ionospheric grid does not match the input grid parameter.");
      return;
    }
    ptr_ionosphere->lon_max = ptr_ionosphere->lon_min + 
		       (ptr_ionosphere->num_lon - 1) * ptr_ionosphere->lon_inc;
    ptr_ionosphere->ht_min = *(mxGetPr(prhs[12]) + 6);
    ptr_ionosphere->ht_inc = *(mxGetPr(prhs[12]) + 7);
    ptr_ionosphere->num_ht = *(mxGetPr(prhs[12]) + 8);
    if (dims[2] != ptr_ionosphere->num_ht) {
      mexErrMsgTxt("The number of heights in the input ionospheric grid does not match the input grid parameter.");
      return;
    }
    ptr_ionosphere->ht_max = ptr_ionosphere->ht_min + 
		      (ptr_ionosphere->num_ht - 1) * ptr_ionosphere->ht_inc;

    dims = mxGetDimensions(prhs[13]);
    for (kk = 0; kk < dims[2]; kk++) {
      for (jj = 0; jj < dims[1]; jj++) {
	for (ii = 0; ii < dims[0]; ii++) {
	  idx = kk*dims[1]*dims[0] + jj*dims[0] + ii;
	  geomag_field.Bx[kk][jj][ii] = *(mxGetPr(prhs[13]) + idx);
	  geomag_field.By[kk][jj][ii] = *(mxGetPr(prhs[14]) + idx);
	  geomag_field.Bz[kk][jj][ii] = *(mxGetPr(prhs[15]) + idx);
	}
      }
    }
    geomag_field.lat_min = *(mxGetPr(prhs[16]) + 0);
    if (geomag_field.lat_min < -90.0 || geomag_field.lat_min > 90.0) {
      mexErrMsgTxt("The start latitude of the input geomagnetic field grid must be in the range -90 to 90 degrees.");
      return;
    }
    geomag_field.lat_inc = *(mxGetPr(prhs[16]) + 1);
    geomag_field.num_lat = *(mxGetPr(prhs[16]) + 2);
    if (dims[0] != geomag_field.num_lat) {
      mexErrMsgTxt("The number of latitudes in the input geomagnetic field grid does not match the input grid parameter.");
      return;
    }
    geomag_field.lat_max = geomag_field.lat_min + 
			   (geomag_field.num_lat - 1) * geomag_field.lat_inc;
    geomag_field.lon_min = *(mxGetPr(prhs[16]) + 3);
    if (geomag_field.lon_min < -180.0 || geomag_field.lon_min > 180.0) {
      mexErrMsgTxt("The start longitude of the input geomagnetic field grid must be in the range -180 to 180 degrees.");
      return;
    }
    geomag_field.lon_inc = *(mxGetPr(prhs[16]) + 4);
    geomag_field.num_lon = *(mxGetPr(prhs[16]) + 5);
    if (dims[1] != geomag_field.num_lon) {
      mexErrMsgTxt("The number of longitudes in the input geomagnetic field grid does not match the input grid parameter.");
      return;
    }
    geomag_field.lon_max = geomag_field.lon_min + 
			   (geomag_field.num_lon - 1) * geomag_field.lon_inc;
    geomag_field.ht_min = *(mxGetPr(prhs[16]) + 6);
    geomag_field.ht_inc = *(mxGetPr(prhs[16]) + 7);
    geomag_field.num_ht = *(mxGetPr(prhs[16]) + 8);
    if (dims[2] != geomag_field.num_ht) {
      mexErrMsgTxt("The number of heights in the input geomagnetic field grid does not match the input grid parameter.");
      return;
    }
    geomag_field.ht_max = geomag_field.ht_min + 
			  (geomag_field.num_ht - 1) * geomag_field.ht_inc;

    /* Now the ionosphere has been read in set the iono_exist_in_mem flag to
       indicate this for future raytrace calls */
    iono_exist_in_mem = 1;
  }
   
  /* If we are doing a magneto-ionic raytrace then check that the magnetic
     field grid is consistent with the ionospheric grid */
  if (OX_mode != 0) {
    if (geomag_field.lat_min != ptr_ionosphere->lat_min) {
      snprintf(error_mesg, 200, "\nThe minimum latitude of the geomagnetic field grid (%f deg.) is inconsistent with the electron density grid (%f deg).", geomag_field.lat_min, ptr_ionosphere->lat_min);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (geomag_field.lat_max != ptr_ionosphere->lat_max) {
      snprintf(error_mesg, 200, "\nThe maximum latitude of the geomagnetic field grid (%f deg.) is inconsistent with the electron density grid (%f deg).", geomag_field.lat_max, ptr_ionosphere->lat_max);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (geomag_field.lon_min != ptr_ionosphere->lon_min) {
      snprintf(error_mesg, 200, "\nThe minimum longitude of the geomagnetic field grid (%f deg.) is inconsistent with the electron density grid (%f deg).", geomag_field.lon_min, ptr_ionosphere->lon_min);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (geomag_field.lon_max != ptr_ionosphere->lon_max) {
      snprintf(error_mesg, 200, "\nThe maximum longitude of the geomagnetic field grid (%f deg.) is inconsistent with the electron density grid (%f deg).", geomag_field.lon_max, ptr_ionosphere->lon_max);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (geomag_field.ht_min  != ptr_ionosphere->ht_min) {
      snprintf(error_mesg, 200, "\nThe minimum height of the geomagnetic field grid (%f deg.) is inconsistent with the electron density grid (%f deg).", geomag_field.ht_min, ptr_ionosphere->ht_min);
      mexErrMsgTxt(error_mesg);
      return;
    }
    if (geomag_field.ht_max  != ptr_ionosphere->ht_max) {
      snprintf(error_mesg, 200, "\nThe maximum height of the geomagnetic field grid (%f deg.) is inconsistent with the electron density grid (%f deg).", geomag_field.ht_max, ptr_ionosphere->ht_max);
      mexErrMsgTxt(error_mesg);
      return;
    }
  }

  /* If a user specified input state vector is not being used then check to make
     sure that start height of raytracing is below the start of the ionosphere.
     If not then quit with error message. */
  if (nrhs != 10 && nrhs != 18) {
    if (start_height > ptr_ionosphere->ht_min) {
      mexErrMsgTxt("The start height for ray tracing must be below the start height of the ionospheric grid. If you want to start ray tracing inside the ionosphere then you must also specify the initial state vector of the ray");
      return;
    }
  } 

  /* read in the optional input (structure containing a user defined starting 
     ray state vector for each ray) from  Matlab (if required) and check to make
     sure each field is valid */
  mxArray *tmp_ptr;
  ray_state_vec_in = (double *) mxMalloc(12*num_rays*sizeof(double));
  if (nrhs == 10 || nrhs == 18) {
    for (kk = 0; kk < num_rays; kk++) {
      for (jj = 0; jj < 12; jj++) {
	idx = kk*12 + jj;
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
	if (mxIsNaN(*(ray_state_vec_in +idx)) ||
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
      for (jj = 0; jj < 12; jj++) {
	idx = kk*12 + jj;
        *(ray_state_vec_in + idx) = -1;
      }
    }
  }
    
  /* init pointers and malloc space for return data */
  nhops_attempted = (int *) mxMalloc(num_rays*sizeof(int));
  npts_in_ray = (int *) mxMalloc(num_rays*sizeof(int));
  ray_data = (double *) mxMalloc(15*nhops*num_rays*sizeof(double));
  ray_path_data =(double *) mxMalloc(19*max_pts_in_ray*num_rays*sizeof(double));
  ray_label = (int *) mxMalloc(nhops*num_rays*sizeof(int));
  ray_state_vec_out =
                (double *) mxMalloc(12*max_pts_in_ray*num_rays*sizeof(double));

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

  /* call the computational routine raytrace_3d_ */
  step_size_min = step_size_min * 1000.0;         /* convert to meters */
  step_size_max = step_size_max * 1000.0;         /* convert to meters */

  raytrace_3d_(&start_lat, &start_long, &start_height, &num_rays, elevs,
	       bearings, freqs, &OX_mode, &nhops, &step_size_min,
	       &step_size_max, &tol, ptr_ionosphere, &geomag_field,
	       ray_state_vec_in, &return_ray_path_data, &return_ray_state_vec,
               ray_data, ray_path_data, ray_label, nhops_attempted, npts_in_ray,
               ray_state_vec_out, &elapsed_time);

  /* Copy the raytrace data into the matlab data structures. */
  /* If only 1 ray has been raytraced (num_rays == 1) and either of the 
     environment variables PHARLAP_OLD_FORMAT or PHARLAP_OLD_FORMAT_3D has been 
     set to "true" then return the data in Matlab Double Matricies. This is to 
     retain backwards compatibility with Matlab code using PHaRLAP version 
     3.7.1 and earlier. */
  new_format = 1;               /* new format is the default */
  old_format = getenv("PHARLAP_OLD_FORMAT");
  old_format_3D = getenv("PHARLAP_OLD_FORMAT_3D");
  
  if (old_format != NULL) {
    if (strcmp(old_format, "true") == 0 && num_rays == 1) new_format = 0;
  }
  if (old_format_3D != NULL) {
    if (strcmp(old_format_3D, "true") == 0 && num_rays == 1) new_format = 0;
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
	idx = kk + num_rays*ray_data_fieldname_raytrace_output_position[jj];
	tmp_ptr = mxCreateDoubleMatrix(1, nhops_attempted[kk], mxREAL);
	mxSetFieldByNumber(plhs[0], (mwIndex) kk, (mwIndex) jj, tmp_ptr);
	stepmemcpyd(mxGetPr(tmp_ptr), &ray_data[idx], num_rays*15,
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
	mxSetField(plhs[1], (mwIndex) kk, "initial_bearing", tmp_ptr);
	*mxGetPr(tmp_ptr) = bearings[kk];

	tmp_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxSetField(plhs[1], (mwIndex) kk, "frequency", tmp_ptr);
	*mxGetPr(tmp_ptr) = freqs[kk];

	for (jj = 0; jj < ray_path_data_numfields-3; jj++) {
	  idx = kk + num_rays*jj;
	  tmp_ptr = mxCreateDoubleMatrix(1, npts_in_ray[kk], mxREAL);
	  mxSetFieldByNumber(plhs[1], (mwIndex) kk, (mwIndex) jj+3, tmp_ptr);
	  stepmemcpyd(mxGetPr(tmp_ptr), &ray_path_data[idx], num_rays*19,
							       npts_in_ray[kk]);
	}
      }

      /* copy the ray_state_vec structure if requested */
      if (return_ray_state_vec) {
        for (jj = 0; jj < ray_state_vec_numfields; jj++) {
	  idx = kk + num_rays*jj;
	  tmp_ptr = mxCreateDoubleMatrix(1, npts_in_ray[kk], mxREAL);
	  mxSetFieldByNumber(plhs[2], (mwIndex) kk, (mwIndex) jj, tmp_ptr);
	  stepmemcpyd(mxGetPr(tmp_ptr), &ray_state_vec_out[idx], num_rays*12,
		                                               npts_in_ray[kk]);
        }
      }

    }
    
  } else {          /* copy data to Matlab Double Matricies */
    
    /* Create matricies for the return arguments to matlab */
    plhs[0] = mxCreateDoubleMatrix(15, nhops_attempted[0], mxREAL);
    plhs[1] = mxCreateDoubleMatrix(19, npts_in_ray[0], mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, nhops_attempted[0], mxREAL);
    plhs[4] = mxCreateDoubleMatrix(12, npts_in_ray[0], mxREAL);

    /* Load the output data from raytrace into the LHS pointers, which are the 
       outputs to MATLAB. */
    for (ii=0; ii< nhops_attempted[0]; ii++) {
      idx = ii*15;
      stepmemcpyd(mxGetPr(plhs[0])+idx, &ray_data[ii], nhops, 15);
    }
    for (ii=0; ii< npts_in_ray[0]; ii++) {
      idx = ii*19;
      stepmemcpyd(mxGetPr(plhs[1])+idx, &ray_path_data[ii], max_pts_in_ray, 19);
    }
    *mxGetPr(plhs[2]) = nhops_attempted[0];

    for (jj = 0; jj < nhops_attempted[0]; jj++) {
      *(mxGetPr(plhs[3]) + jj) = (double) (&ray_label[0])[jj];
    }

    if (nlhs == 5) {
      for (ii=0; ii< npts_in_ray[0]; ii++) {
        idx = ii*12;
        stepmemcpyd(mxGetPr(plhs[4])+idx, &ray_state_vec_out[ii],
		                                        max_pts_in_ray, 12);
      }
    }
    
  }

  /* free malloced space */
  mxFree(elevs);
  mxFree(bearings);
  mxFree(freqs);
  mxFree(ray_state_vec_in);
  mxFree(nhops_attempted);
  mxFree(npts_in_ray);
  mxFree(ray_data);
  mxFree(ray_path_data);
  mxFree(ray_label);
  mxFree(ray_state_vec_out);

}



