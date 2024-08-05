%
% Name:
%   ART_QP.m
%
% Purpose: 
%   Analytical raytrace for an ionosphere constructed from quasi-parabolic 
%   layers
%
% Calling Sequence:
%   [ground, group, phase] = ART_QP(elevs, freq, re, QPcoeffs);
%
% Inputs:
%   elev     -   elevation of the ray (degrees)
%   freq     -   frequency of the ray (MHz)
%   re       -   radius of the Earth (approximated as spherical) (km)
%   QPcoeffs -   array of QP coefficients defining each QP segment, the
%                routine QP_profile.m generates this array
%     QPcoeffs(ii, 1) = a     where a = Nm 
%     QPcoeffs(ii, 2) = b     where b = Nm .* ((hm + re - ym) / ym).^2
%     QPcoeffs(ii, 3) = rb    radial coordinate of bottom of segment
%     QPcoeffs(ii, 4) = rt    radial coordinate of top of segment
%     QPcoeffs(ii, 5) = 1 for ionospheric layer (E, F1 or F2)
%                      -1 for joining layer (bottom to E, E to F1, or F1 to F2)
%
% Outputs:
%   ground  -   ground range of landing point ray from launch point (km)
%               NaN returned if ray penetrates the ionosphere
%   group   -   group range of landing point ray from launch point (km)
%               NaN returned if ray penetrates the ionosphere
%   phase   -   phase path of landing point ray from launch point (km)
%               NaN returned if ray penetrates the ionosphere
%
% Dependencies:
%   None
%
% References:
%   Croft and Hoogasian, Radio Science, Vol. 3, 69--74, 1968
%   Dyson and Bennett, J. Atmos. Terr. Phys., Vol. 50, 251--262, 1988
%   Bennett, Chen and Dyson, Applied Computational Electromagnetic Society 
%                            Journal, Vol. 56, 192-210, 1991
%   
% Modification History
%   M. A. Cervera  30/08/2016  V1.0
%      Initial version.
%

function [ground, group, phase] = ART_QP(elev, freq, re, QP_seg_coeffs)

% constants
dtor = pi ./ 180;
plas_fac = 80.6163849431291d-6;   % converts plasma freq to electrons / cm^3
 
% other parameters
rb = QP_seg_coeffs(1, 3);      % bottom of first segment == bottom of ionosphere
elev_r = elev .* dtor;
cos_elev = cos(elev_r);
sin_elev = sin(elev_r);
cos_gamma = re .* cos_elev ./ rb;
sin_gamma = sqrt(1 - cos_gamma.^2);
gamma_r = acos(cos_gamma);
N = freq.^2 ./ plas_fac;       % convert operating frequency to an equivalent
                               % electron density (electrons / cm^3)
sz = size(QP_seg_coeffs);
num_segments = sz(1);

% loop over the QP segments and calculate the Dyson and Bennett integrals
ground_integral = 0.0;
group_integral = 0.0;
phase_integral = 0.0;
apogee_segment = 0;
penetrated_iono = 0;
seg_idx = 0;
while (~apogee_segment & ~penetrated_iono)
  seg_idx = seg_idx + 1;
  
  % QP segment coefficents
  a = QP_seg_coeffs(seg_idx, 1);
  b = QP_seg_coeffs(seg_idx, 2);
  rb = QP_seg_coeffs(seg_idx, 3);           % bottom of segment
  rt = QP_seg_coeffs(seg_idx, 4);           % top of segment
  
  
  % Calculate the Croft and Hoogasian A, B, C coefficients (see page 72 of the
  % reference)
  if QP_seg_coeffs(seg_idx, 5) == 1
    rc = rt;       % this QP segment is an ionospheric layer (E, F1 or F2)
  else
    rc = rb;       % this reverse QP segment joins ionospheric layers 
  end  
  A = 1 - a./N + b./N;                 
  B = -2 .* rc .* b./N;              
  C = (rc.^2 .* b)./N - re.^2 .* cos_elev.^2;  
  Bsq_minus_4AC = B.^2 - 4.*A.*C;      
  
  % determine if this segment is the final reflecting segment
  if (Bsq_minus_4AC > 0) 
    r_apogee = (-B - sqrt(Bsq_minus_4AC))./(2.*A);
    if r_apogee < rt
      apogee_segment = 1;    % ray is reflected in this segment
    else
      apogee_segment = 0;    % ray penetrates this segment
    end
  else 
    apogee_segment = 0;      % ray penetrates this segment
  end
  
  % determine if the ionosphere has been penetrated
  if (seg_idx == num_segments & apogee_segment == 0)
    penetrated_iono = 1;  % ray penetrated final segment and thus ionosphere
  else
    penetrated_iono = 0;
  end

  % Calculate the Dyson and Bennett integrals for ground range - see Dyson
  % and Bennett equation 12 and appendix A4 equations A19 - A22. See also
  % Croft and Hoogasian equation 6a.
  R_L = A.*rb.^2 + B.*rb + C;
  R_U = A.*rt.^2 + B.*rt + C;
  if C > 0
    if apogee_segment
      % Dyson and Bennett equation A21 for top of final reflecting layer
      I1_U = -log(Bsq_minus_4AC) ./ (2.*sqrt(C));
    else
      % Dyson and Bennett equation A19 for top of layer
      I1_U = -log(abs(2.*sqrt(C.*R_U) + B.*rt + 2.*C) ./ rt) ./ sqrt(C);    
    end
    % Dyson and Bennett equation A19 for bottom of layer
    I1_L = -log(abs(2.*sqrt(C.*R_L) + B.*rb + 2.*C) ./ rb) ./ sqrt(C);
  else
    if apogee_segment
      % Dyson and Bennett equation A22 for top of final reflecting layer
      I1_U = pi ./ (2.*sqrt(-C));      
    else
      % Dyson and Bennett equation A20 for top of layer
      I1_U = asin((B.*rt + 2.*C)./(abs(rt).*sqrt(Bsq_minus_4AC))) ./ sqrt(-C);
    end
    % Dyson and Bennett equation A20 for bottom of layer
    I1_L = asin((B.*rb + 2.*C)./(abs(rb).*sqrt(Bsq_minus_4AC))) ./ sqrt(-C);
  end
  
  ground_integral = ground_integral + I1_U - I1_L;
  
  % Calculate the Dyson and Bennett integrals for group range - see Dyson and
  % Bennett eqn 13 and and Appendix A4 equations A23 - A25. See also
  % Croft and Hoogasian equation 6b. Calculate  the Dyson and Bennett
  % integrals for phase path - see Bennett, Chen and Dyson 1991.
  if A > 0
    if apogee_segment
      % Dyson and Bennett equation A24 (A>0) for top of layer
      I3_U = log(Bsq_minus_4AC) ./ (2.*sqrt(A));
    else
      % Dyson and Bennett equation A23 (A>0) for top of layer
      I3_U = log(abs(2.*sqrt(A.*R_U) + 2.*A.*rt + B)) ./ sqrt(A);
    end
    % Dyson and Bennett equation A23 (A>0) for bottom of layer
    I3_L = log(abs(2.*sqrt(A.*R_L) + 2.*A.*rb + B)) ./ sqrt(A);
  else
    if apogee_segment
      % Dyson and Bennett equation A25 (A<0) for top of layer
      I3_U = pi ./ (2.*sqrt(-A));
    else
      % Dyson and Bennett equation A23 (A<0) for top of layer
      I3_U = -asin((2.*A.*rt + B)./ sqrt(Bsq_minus_4AC)) ./ sqrt(-A);
    end
    % Dyson and Bennett equation A23 (A<0) for bottom of layer
    I3_L = -asin((2.*A.*rb + B)./ sqrt(Bsq_minus_4AC)) ./ sqrt(-A);
  end  
  if apogee_segment
    % Dyson and Bennett equation A24 for top of layer
    I2_U = -B.*I3_U ./ (2.*A);
  else
    % Dyson and Bennett equation A23 for top of layer
    I2_U = sqrt(R_U)./A - B.*I3_U ./ (2.*A);
  end
  % Dyson and Bennett equation A23 for bottom of layer
  I2_L = sqrt(R_L)./A - B.*I3_L ./ (2.*A);
  
  group_integral = group_integral + I2_U - I2_L;
  phase_integral = phase_integral + (I3_U - I3_L).*B + (I2_U - I2_L).*A + ...
                                    (I1_U - I1_L).*(C + (re.*cos_elev).^2);
  	
end  % of while (~apogee_segment & ~penetrated_iono)

% Calculate ground and group ranges - equations 12 and 13 of Dyson and Bennett
% Return NaNs if ray has penetrated the final layer.
if penetrated_iono
  ground = NaN;
  group = NaN;
  phase = NaN;
else   
  rb = QP_seg_coeffs(1, 3);  % bottom of first segment == bottom of ionsophere
  ground = 2.*re.*(gamma_r - elev_r + re.*cos_elev.*ground_integral);
  group = 2.*(rb.*sin_gamma - re.*sin_elev + group_integral);
  phase = 2.*(rb.*sin_gamma - re.*sin_elev + phase_integral);
end

end
