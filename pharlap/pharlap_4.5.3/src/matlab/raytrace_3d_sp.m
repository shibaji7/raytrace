% 
% Name : 
%   raytrace_3d_sp
%
% Purpose :
%   3D magneto-ionic numerical raytrace for multi-hop rays. The computational
%   routine is multi-threaded and accepts a user specified number of independant
%   rays which it farms out to the available computational cores.
%
%   The Hamiltonian for the ray is that set down by Haselgrove and Hasegrove
%   (1960), and Hasegrove (1963). Geomagnetic field effects are considered
%   for O and X  polarized radio waves. No-field case is also considered if
%   required.  Spherical Earth coordinate system is assumed with user input
%   radius of the Earth. This version takes approximately 70% of the time to
%   complete that the WGS84 compliant raytrace_3d does. An appropriate choice
%   for   the radius of the Earth will minimise the error wrt the WGS84
%   coordinate  system (see also earth_radius_wgs84.m)
%
% Calling sequence :
%   The first time raytrace_3d_sp is called, ionospheric and geomagnetic grids
%   must be specified. For subsequent calls if these grids are not passed in 
%   (see calling sequence examples 4 5 and 6) then the grids from the
%   previous call (which are in memory) will be used. This may speed up the 
%   subsequent calls to raytrace for large grids. 
%
%   1.  [ray_data, ray_path_data, ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elevs, ray_bearings,
%                       freqs, OX_mode, nhops, tol, rad_earth, iono_en_grid, ...
%    	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
%                       Bx, By, Bz, geomag_grid_parms)
%
%   2.  [ray_data, ray_path_data, ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elevs,ray_bearings,
%                       freqs, OX_mode, nhops, tol, rad_earth, iono_en_grid, ...
%    	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
%                       Bx, By, Bz, geomag_grid_parms, ray_state_vec_in)
%
%   3.  [ray_data, ray_path_data, ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elevs, ray_bearings,
%                       freqs, OX_mode, nhops, tol, rad_earth)
%
%   4.  [ray_data, ray_path_data, ray_state_vec, NRT_execution_time] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elevs, ray_bearings,
%                       freqs, OX_mode, nhops, tol, rad_earth, ray_state_vec_in)
%
% Inputs :
%   origin_lat      - geocentric latitude (-90 to 90 degrees) of start point
%                     of rays (All rays have the same origin)
%   origin_long     - geocentric longitude (-180 to 180 degrees) of start
%                     point of rays
%   origin_height   - height above sea-level of start point of ray which
%                     must be below start of ionosphere (km). (If the start of 
%                     the ray tracing is inside the ionosphere then the initial
%                     state vector of the ray must be specified.)
%   elevs           - 1 X M (where M is the number of rays) array of initial
%                     elevation of rays (deg) 
%   ray_bearings    - 1 X M array of initial bearing of the rays (deg) 
%   freqs           - 1 X M array of wave frequency of the rays (MHz)
%   OX_mode         - polarization mode of ray: 1 = O, -1 = X, 0 = no field
%   nhops           - number of hops to complete   
%   tol             - a 1 or 3 element vector controlling ODE solver precision.
%      tol(1)  =  ODE solver tolerence, valid values 1e-12 to 1e-2
%      tol(2)  =  ODE solver minimum step size to consider (0.001 to 1 km)
%      tol(3)  =  ODE solver maximum step size to consider (1 to 100 km)
%      
%      If tol is a scalar then min and max step sizes are set to 0.01 and
%      10km. Suggested value for tol in this case is 1e-7. This provides
%      backward compatibility with PHaRLAP version 3.2.1 and earlier.
%
%      Example values: 1. Highest precision, slowest  - tol(1) = 1e-8
%                                                       tol(2) = 0.01 km
%                                                       tol(2) = 10 km
%                      2. Lower precision, faster     - tol(1) = 1e-7
%                                                       tol(2) = 0.025 km
%                                                       tol(2) = 25 km 
%                      3. Lowest precision, fastest   - tol(1) = 1e-6
%                                                       tol(2) = 0.1 km
%                                                       tol(2) = 100 km 
%    
%      The example values may be selected by setting tol to be a scalar with
%      a value of 1, 2, or 3.
%
%   rad_earth       - radius of Earth to use for coord conv. routines (m)
%   iono_en_grid    - 3d grid of ionospheric electron density (electrons/cm^3).
%                     Maximum grid size is 701 X 701 X 401 (lat, lon, height).
%                     The ionospheric data and their gradients must be smooth
%                     and continuous. See also note 2 below on memory.  
%   iono_en_grid_5  - 3d grid of ionospheric electron density (electrons/cm^3)
%                     5 minutes later. This array is used to calculate Doppler
%                     shift. If Doppler shift is not required then this array
%                     can be set to iono_en_grid. Size of grid must be
%                     identical to iono_en_grid. 
%   collision_freq  - 3d grid of effective electron collision frequency (Hz). 
%                     This is used to calculate the total and deviative
%                     absorption (see notes). If absorption is not required
%                     then the array elements can be set to zero. Size of
%                     grid must be identical to iono_en_grid. 
%
%   iono_grid_parms - 9x1 vector containing the parameters which define the
%                     ionospheric grid :
%           (1) geocentric latitude (degrees) of start of grid - must be in the
%               range -90 to 90 degrees 
%           (2) latitude step (degrees)
%           (3) number of latitudes
%           (4) geocentric longitude (degrees) of start of grid - must be in
%               the range -180 to 180 degrees
%           (5) lonfitude step (degrees)
%           (6) number of longitudes
%           (7) geocentric height (km) start
%           (8) height step (km)
%           (9) number of heights
%
%   Bx, By, Bz        - 3D grids of x, y and z components of the geomagnetic 
%                       field (Tesla) as a function of latitude,
%                       longitude and height. The maximum grid size is 
%                       101 X 101 X 201 (lat, lon, height). The geomagetic field
%                       data and their gradients must be smooth and continuous.
% 
%   geomag_grid_parms - 9x1 vector containing the parameters which define the
%                       ionospheric grid :
%           (1) geocentric latitude (degrees) of start of grid - must be in the
%               range -90 to 90 degrees 
%           (2) latitude step (degrees)
%           (3) number of latitudes
%           (4) geocentric longitude (degrees) of start of grid - must be in
%               the range -180 to 180 degrees
%           (5) longitude step (degrees)
%           (6) number of longitudes
%           (7) geocentric height (km) start
%           (8) height step (km)
%           (9) number of heights
%
% Optional Inputs :
%   ray_state_vec_in(:) - 1 X M structure containing the state vector of each
%            ray (total of M rays) at its starting point. Each field (total of
%            12 fields) is an element of the state vector. If the first field
%            for a given ray is -1 then default starting values for a ray
%            launched from the ground are used. NB: non-default values are
%            for EXPERTS ONLY - you really want to know what you are doing!
%            If used then the input starting point, elevation and bearing of
%            ray (first 5 inputs) will be ignored. The fields of this structure
%            are:
%     .pos_x          - cartesian x position of ray start point (m)
%     .pos_y          - cartesian y position of ray start point (m)
%     .pos_z          - cartesian z position of ray start point (m)
%     .dir_x          - cartesian x, y, and z components of vector with a
%     .dir_y            magnitude equal to the refractive index pointing in the
%     .dir_z            direction of the wave normal (which, generally, is not 
%                       in the direction  of the ray propagation) 
%     .group_path     - group path (m)
%     .geom_path      - geometrical distance travelled by rays (m)
%     .phase_path     - phase path (m)
%     .absorption     - ionospheric absorption (dB)
%     .indep_var      - independant variable 
%     .ODE_step_size  - ODE solver independant variable step size
%                                                                            
% Outputs :
%   ray_data(:) -  1 X M structure containing the information for each of the
%                  rays (total of M rays). Each field is a 1 X N array
%                  containing information for each hop. The fields of this
%                  structure are: 
%     .lat                   - geocentric latitude of end point of ray (deg) 
%                              for each hop 
%     .lon                   - geocentric longitude of end point of ray (deg)
%                              for each hop 
%     .ground_range          - geocentric ground range (Km) 
%     .group_range           - group range (Km) for each hop    
%     .phase_path            - phase path (km) for each hop
%     .initial_elev          - initial elevation (deg) at this hop   
%     .final_elev            - final elevation (deg) of this hop, IEEE NaN
%                              value is assigned if ray does not return to
%                              ground         
%     .initial_bearing       - initial ray bearing (deg) at this hop 
%     .final_bearing         - final ray bearing (deg) of this hop 
%     .total_absorption      - total ionospheric absorption (dB) (see notes)   
%                              cumulative over hops
%     .deviative_absorption  - ionospheric deviative absorption (dB) (see notes)   
%                              cumulative over hops
%     .TEC_path              - total electron content along ray path (number
%                              of electrons in 1m^2 cross-section tube)  
%     .doppler_shift         - Doppler shift (Hz)
%     .apogee                - apogee of the ray (km)
%     .geometric_path_length - Geometrical distance travelled by ray (km)
%     .frequency             - carrier frequency of the ray (MHz)
%     .nhops_attempted       - number of hops actually attempted
%     .NRT_execution_time    - time taken for computer to complete raytracing (s)
%     .ray_label             - label for each hop attempted which indicates
%                              what the ray has done. 
%          = 1    for ray reaching ground                           
%           -1    for field aligned backscatter - ray reflected with appropriate
%                 scattering loss, raytracing terminated
%                 *** NB FAI BACKSCATTER NOT YET IMPLEMENTED ***
%           -2    ray has penetrated the ionosphere - raytracing terminated 
%           -3    ray has exited ionospheric grid - raytracing terminated
%           -4    ray has exceeded the maximum allowed points along path
%                 (20000 points) raytracing terminated
%           -100  a catastrophic error occured - terminate raytracing
%
%   ray_path_data(:) - 1 X M structure containing information about each of the
%                      rays at each point along their paths. The cartesian
%                      coordinate system (below) is defined with the x axis
%                      passing through the equator at the prime meridian, y
%                      through the equator at longitude of +90 degrees, and z
%                      through the geographic north pole.
%     .initial_elev            - initial elevation of the ray (degrees)
%     .initial_bearing         - initial bearing of the ray (degrees)
%     .frequency               - carrier frequency of the ray (MHz)
%     .lat                     - geocentric latitude (degrees)
%     .lon                     - geocentric longitude (degrees)
%     .height                  - height of ray (km)
%     .group_range             - group path (km)
%     .phase_path              - phase path (km)
%     .refractive_index        - refractive index
%     .group_refractive_index  - group refractive index
%     .wavenorm_ray_angle      - angle between wave normal and the ray 
%                                direction (degrees) 
%     .wavenorm_B_angle        - angle between wave normal and geomagnetic
%                                field (degrees)
%     .polariz_mag             - magnitude of the wave volume-polarization
%                                vector, R (see Notes)
%     .wave_Efield_tilt        - wave E-field tilt angle out of the plane of 
%                                the wave front 
%     .volume_polariz_tilt     - volume polarization vector tilt angle out of
%                                the wave-front plane. NB. If e- = 0 then the 
%                                magnitude of the polarization field is zero. In
%                                this case having a polarization field tilt does
%                                not make sense and so it is set IEEE NaN.
%     .electron_density        - electron density (e- / cm^3)
%     .geomag_x                - WGS84 x component of geomagnetic field (Tesla)
%     .geomag_y                - WGS84 y component of geomagnetic field (Tesla)
%     .geomag_z                - WGS84 z component of geomagnetic field (Tesla)
%     .geometric_distance      - geometrical distance travelled by ray (km)
%     .collision_frequency     - collision frequency at each point along ray
%     .cumulative_absorption   - cumulative_absorption along ray path (dB)
%                                                                           
%   ray_state_vec(:) - 1 X M structure containing the state vector of each
%                      ray at each point along their paths. The cartesian
%                      coordinate system (below) is defined with the x axis
%                      passing through the equator at the prime meridian, y
%                      through the equator at longitude of +90 degrees, and z
%                      through the geographic north pole.
%      .pos_x            - cartesian x position of ray (m)
%      .pos_y            - cartesian y position of ray (m)
%      .pos_z            - cartesian z position of ray (m)
%      .dir_x            - cartesian x, y, and z components of vector 
%      .dir_y              with a magnitude equal to the refractive index  
%      .dir_z              pointing in the direction of the wave normal
%      .group_path       - group path (m)
%      .geometrical_path - geometrical distance travelled by rays (m)
%      .phase_path       - phase path (m)
%      .absorption       - ionospheric absorption (dB)
%      .indep_var        - independant variable 
%      .ODE_step_size    - ODE solver independant variable step size
%
%
% Notes:
% 1. Total and Deviative absorption
%   The total absorption is calculated by integrating the imaginary component
%   of the complex refractive index along the ray path. See: 
%   Pederick and Cervera (2014), Radio Sci., 49, 81-93,doi:10.1002/2013RS005274 
%   Here we use the Appleton-Hartree formulation of the complex refractive
%   index which requires the effective collision frequency to be specified in
%   the collision frequency grid. See:
%   Zawdie et al. (2017), Radio Sci., 52, 767-783, doi:10.1002/2017RS006256.
%
%   Deviative absorption occurs where ray bending or group retardation takes
%   place. Here we calculate the deviative absorption by integrating equation
%   (9) of Appleton and Piggott (1954) above the E-layer (above an altitude
%   of hmE + ymE). See e.g. :
%   Appleton and Piggott (1954), J. Atmos. Terr. Phys., 5, 141-172, 
%      doi:10.1016/0021-9169(54)90029-X
%   Davies (1990), Ionospheric Radio, pp214-261
%   Coleman (1997), J. Atmos. Solar-Terrestrial. Phys., 59, 2089-2099,
%      doi:10.1016/S1364-6826(97)00038-2
%
%   Deviative absorption will only be required to be used if an alternate
%   technique is used to calculate absorption which does not adequately
%   account for it. George and Bradley absorption accounts for deviative
%   absorption  in the E layer but not in the F layer. Thus, if abso_bg.m is
%   used to calculate absorption, deviative absorption will need to be added.
%
% 2. Polarization :
%   The formulae to calculate the wave polarization are those from Davies (1990)
%   pp78 and assume that collisions are negligible, for reference see also 
%   Budden (1985), pp66-74. The coordinate system (cartesian i, j, k) is defined
%   as follows : the wave-normal is in the i-direction and the 
%   geomagnetic field direction lies in the i-j plane (thus k is orthoganal 
%   to B-field). Then the magnitude of the wave polarization, R, gives the
%   axial ratio of the ellipse in the j-k plane that is traced out by the 
%   electric field vector of the wave projected onto the j-k plane. The 
%   electric field vector of the wave actually lies in the plane tilted  
%   forward (or backward) by angle polariz_E_ang from the j axis for the  
%   O-mode (or X-mode). For an O-mode (or X-mode) ray, the semi-major axis of
%   the polarization ellipse projected in the j-k plane is in the j-direction
%   (or k-direction). See Figures 3.1 and 3.3 of Davies (1990).
%
% 3. Memory usage :
%   If large ionospheric grid sizes are used then the required memory could
%   exceed that available. In this case memory paging will occur which will
%   greatly slow down the computation speed. If the maximum allowable size
%   for the ionospheric grids is specified, then ~4.4GB of memory is required
%   (~1.8GB on 32-bit Windows) in addition to all the other memory
%   requirements (701 X 701 X 401 elements X 3 iono grids X 8 bytes). 
%
% 4. Multi-threading :
%   Raytrace calls with mulitple rays will automatically have the separate
%   rays distributed over several computational threads. The number of
%   threads defaults to the number of processor cores available. Setting the
%   environement variable NRT_NUM_THREADS will override this - e.g.
%       >> setenv('NRT_NUM_THREADS', 6) 
%   will limit the number of computational thread to 6 even if more cpu cores
%   are available. Valid values are 1 to 64. Invalid values will result in
%   the default being used.
%
% References:
%   1. Haselgrove, C.B., and Haselgrove, J. (1960), "Twisted ray paths in the
%      ionosphere", Proc. Phys. Soc. (London), Vol. 75, 357-361.
%
%   2. Haselgrove, J. (1963), "The Hamiltonian ray path equations", J. Atmos. 
%      Terr. Phys., Vol. 25, 397-399. 
%  
%   3. Davies, K. (1990), "Ionospheric Radio", IEE Electromagnetic Waves 
%      Series 31, Peter Peregrinus, London.
%
%   4. Budden, K. G. (1985), "The propagation of radio waves", Cambridge 
%      University Press, Cambridge.
%
%   5. Bennett, J. A. (1967), "The calculation of Doppler Shifts due to a 
%      changing ionosphere", J. Atmos. Terr. Phys., 1967, Vol. 29, 887-891.
%

% This a Matlab help file only. The actual programme is a mex wrapper
% (raytrace-3d_matlab_wrapper.for) to the Fortran code (raytrace_3d.for).
%
% Modification history:
%   19-06-2010 M.A.Cervera  Initial version.
%                           
%   02-12-2015 M.A.Cervera
%      Updated for multiple ray input/ouput capabilities of PHaRLAP 4.0.0
%                           
%   04-05-2017 M. A. Cervera
%     Implemented a more efficient algorithm for calculating the partial 
%     derivatives of the electron density and the components of the geomagnetic
%     field
%
