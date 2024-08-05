%
% Name:
%   QP_profile_sl.m
%
% Purpose:
%   Returns electron density profile for a single QP layer with a small
%   reverse QP layer at the bottom of the profile to give a smooth increasing
%   electron density from zero.
%   
%
% Calling Sequence:
%   plasma = QP_profile_sl(fc, hm, ym, heights, re) 
%
% Inputs:
%   fm        -   maximum plasma frequency of the QP layer (MHz)
%   hm        -   height of maximum electron density of QP layer (km)
%   ym        -   semi-thickness of the QP layer
%   height    -   array of heights at which QP profile is to be calculated (km)
%   re        -   radius of the Earth (approximated as spherical) (km)
%
% Outputs:
%   elec_dens     -  electron density profile (elecrons / cm^3)
%   QP_seg_coeffs -  the coefficients of the quasi parabolic segments (both
%                    actual layers and the joining layers) which define the
%                    electron density profile
%
% Dependencies:
%   None
%
% References:
%   Dyson and Bennett, J. Atmos. Terr. Phys., Vol. 50, 251-262, 1988
%   Croft and Hoogasian, Radio Science, Vol. 3, 69--74, 1968
%   
% Modification History
%   M. A. Cervera  13/07/2015  V1.0
%      Initial version.
%

function [elec_dens, QP_seg_coeffs] = QP_profile_sl(foE, hmE, ymE, heights, re) 
 
  plas_fac =  80.6163849431291d-6;   % converts plasma freq to electrons / cm^3
  
  r = heights + re;
    
  % calculate the electron density profile for the E layer
  rm = hmE + re;
  rb = rm - ymE;  
  Nm = foE.^2 ./ plas_fac;
  a = Nm;
  b = Nm .* (rb / ymE).^2;
  N_layer = a - b.*(1 - rm./r).^2;
  idx = find(N_layer < 0);
  N_layer(idx) = 0;
    
  % calculate the inverse QP electron density profile for the joining segment
  % which smooths the bottom of E layer
  aj = 0.0;
  rj = rm - 1.2*ymE;
  rc = (rm .* b .* (rm./rj - 1)) ./ (a - aj + b.*(rm./rj - 1));
  bj = -rm .* b .* (1 - rm./rc) ./ (rj .* (1 - rj./rc));
  N_join = aj + bj.*(1 - rj./r).^2;
  idx = find(heights+re < rj );
  N_join(idx) = 0;
  
  % populate the QP segment coefficeints array for the E layer segment and
  % the smoothing segment
  segment = 2;
  QP_seg_coeffs(segment, 1) = a;
  QP_seg_coeffs(segment, 2) = b;
  QP_seg_coeffs(segment, 3) = rc;
  QP_seg_coeffs(segment, 4) = rm;
  QP_seg_coeffs(segment, 5) = 1;   % E layer segment

  segment = 1;
  QP_seg_coeffs(segment, 1) = aj;
  QP_seg_coeffs(segment, 2) = -bj;
  QP_seg_coeffs(segment, 3) = rj;
  QP_seg_coeffs(segment, 4) = rc;
  QP_seg_coeffs(segment, 5) = -1;   % joining segment - smooths bottom of E


  % now calculate the plasma frequency profile
  plasma = sqrt(N_layer .* plas_fac);
  plasj = sqrt(N_join .* plas_fac);
  
  N_tot = N_layer;
  idx = find(heights+re < rc);
  N_tot(idx) = N_join(idx);
  
  elec_dens = N_tot;
  
end
