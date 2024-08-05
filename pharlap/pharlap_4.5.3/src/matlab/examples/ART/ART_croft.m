%
% Name:
%   ART_croft.m
%
% Purpose: 
%   Analytical Raytrace for a single QP layer - Croft and Hoogasian, Radio Sci,
%   Vol. 3, 69--74, 1968
%
% Calling Sequence:
%   [ground, group, phase] = ART_croft(elevs, freq, fc, hm, ym, re);
%
% Inputs:
%   elev  -   elevation of the ray (degrees)
%   freq  -   frequency of the ray (MHz)
%   fc    -   critical frequency of the QP layer (MHz)
%   hm    -   height of maximum electron density of QP layer (km)
%   ym    -   semi-thickness of the QP layer
%   re    -   radius of the Earth (approximated as spherical) (km)
%
% Outputs:
%   ground  -   ground range of landing point ray from launch point (km)
%               NaN returned if ray penetrates the layer
%   group   -   group range of landing point ray from launch point (km)
%               NaN returned if ray penetrates the layer
%   phase   -   phase path of landing point ray from launch point (km)
%               NaN returned if ray penetrates the layer
%
% Dependencies:
%   None
%
% References:
%   Croft and Hoogasian, Radio Science, Vol. 3, 69--74, 1968
%   
% Modification History
%   M. A. Cervera  06/10/2014  V1.0
%      Initial version.
%

function [ground, group, phase] = ART_croft(elev, freq, fc, hm, ym, re)

% constants
dtor = pi ./ 180;

% other parameters
F = freq / fc;        % ratio of operating freq. to critical freq. of QP layer
rb = re + hm - ym;    % bottom of ionosphere
rm = re + hm;
cos_elev = cos(elev .* dtor);
sin_elev = sqrt(1 - cos_elev.^2);
cos_gamma = re .* cos_elev ./ rb;
sin_gamma = sqrt(1 - cos_gamma.^2);
elev_r = elev .* dtor;
gamma_r = acos(cos_gamma);

A = 1 - 1./F.^2 + (rb ./ (F .* ym)).^2;               % Croft and Hoogasian pp72
B = -2 .* rm .* rb.^2 ./ (F.^2 .* ym.^2);             % Croft and Hoogasian pp72
C = (rb.*rm ./ (F .* ym)).^2 - re.^2 .* cos_elev.^2;  % Croft and Hoogasian pp72
		
sqrt_A = sqrt(A);
sqrt_C = sqrt(C);
Bsq_minus_4AC = B.^2 - 4.*A.*C;

% Calculate ground range (Croft and Hoogasian eqn 6a)
ground = 2.*re .* ( (gamma_r - elev_r) - 0.5.*re.*cos_elev./sqrt_C .* ...
  log(Bsq_minus_4AC ./ (4.*C .* (sin_gamma + sqrt_C./rb + 0.5.*B./sqrt_C).^2)));

% Calculate group range (Croft and Hoogasian eqn 6b)
group = 2 .* (rb.*sin_gamma - re.*sin_elev + (1./A) .* ...
       ( -rb.*sin_gamma - B./(4.*sqrt_A) .* ...
         log(Bsq_minus_4AC ./ (2.*A.*rb + B + 2.*rb.*sqrt_A.*sin_gamma).^2) ...
       ));
   
% calculate phase path (Croft and Hoogasian eqn 6c)
phase = -2.*re*sin_elev + 0.5.*B  .* ...
    ( log(Bsq_minus_4AC ./ ...
                   (4.*(A.*rb + B./2 + sqrt_A.*rb.*sin_gamma).^2) ...
	 ) ./ sqrt_A + ...
      log(Bsq_minus_4AC ./ ...
                   (4.*C.*(sin_gamma + sqrt_C./rb + B./2./sqrt_C).^2) ...
	 ) .* rm ./ sqrt_C ...
    );

% check to see if ray has penetrated
if Bsq_minus_4AC < 0
  ground = NaN;
  group = NaN;
  phase = NaN;
end

end
