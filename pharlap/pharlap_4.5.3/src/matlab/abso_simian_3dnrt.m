%
% Name :
%   abso_simian_3dnrt.m
% 
% Purpose / description:
%   Calculates absorption loss for a given ray from a 3D ray trace, using
%   the SiMIAN absorption model. Note - this will typically be not needed as 
%   the 3D NRT routines raytrace_3d.mex and raytrace_3d_sp.mex return the
%   absorption in their ray_data output structures (see the raytrace_3d.m
%   help file). However, abso_simian_3dnrt.m allows the calculation of the
%   absorption using the Sen-Wyller formulation of the complex refractive
%   index in addition to the Appleton-Hartree formulation. The
%   raytrace_3d.mex and raytrace_3d_sp.mex routines calculate the 
%   absorption using Appleton-Hartree only for speed considerations. The
%   absorption calculation using Sen-Wyller is signaficantly slower. It is
%   noted that the Sen-Wyller formulation of the complex refractive gives
%   a value for the absorption which is very close to that when Appleton-
%   Hartree is used if the appropriate form of the electon collision
%   frequnecy is also used (effective collision frequency for A-H
%   vs. mono-energetic collision frequency for S-W). See Zawdie et al. (2017),
%   Radio Sci, vol. 52, 767-783, doi:10.1002/2017RS006256
%
%   The example routine abso_comp_3dnrt.m compares absorption calculated
%   using different methods including Sen-Wyller.
%
% Calling sequence:
%   absorption = abso_simian_3dnrt(ray_path_data, UT, OX_mode, Sen_Wyller)
%
% Inputs (scalar unless indicated):
%   ray_path_data - array containing information about the ray at each
%                   point along its path, as output by raytrace_3d
%   UT            - 5x1 vector containing UTC date and time - year, month,
%                   day, hour, minute
%   OX_mode       - indicates polarisation mode
%                     1 = O-mode
%                    -1 = X-mode
%                     0 = ignore geomagnetic field
%
%   Sen_Wyller    - flag to indicate whether or not to use the Sen_Wyller
%                   formulation
%                     logical true  = use Sen_Wyller
%                     logical false = use Appleton-Hartree
%
% Outputs:
%   absorption   - total absorption loss for the path (dB)
%
% Author:
%   V1.0  L.H. Pederick 07/06/2012
%       Initial version
%
%   V2.0  M. A. Cervera 17/09/2019
%      Re-based the code to work with PHaRLAP 4.3.0
%

function absorption = abso_simian_3dnrt(ray_path_data, UT, OX_mode, Sen_Wyller)

  speed_light = 2.99792458e8;
  neper_to_dB = 8.685889638065037;    % 20*log10(e)

  % get various quantities along the ray path  
  nu = ray_path_data.collision_frequency(1:end);
  %ds = (ray_path_data.geometric_distance(2:end) - ...
  %	  ray_path_data.geometric_distance(1:end-1)) *1e3; % m
  %ds = [0 ds];
  ds = deriv(ray_path_data.geometric_distance)*1e3;
  N = ray_path_data.electron_density(1:end)*1e6;  % electrons per m^3

  mag_field = sqrt(ray_path_data.geomag_x.^2 + ray_path_data.geomag_y.^2 + ...
		   ray_path_data.geomag_z.^2);
  mag_field = mag_field(1:end) ;
  theta = ray_path_data.wavenorm_B_angle(1:end);
  freq = ray_path_data.frequency * 1e6; 
  ang_freq = 2*pi * ray_path_data.frequency * 1e6; 

  if (OX_mode == 0) mag_field = 0; end

  % Calculate the imaginary component of the complex refractive index 
  % using Sen-Wyller or Appleton-Hartree
  if (Sen_Wyller)  % Sen-Wyller

    % convert effective collision frequency to be suitable to use with 
    % Sen-Wyller
    nu = (2/5)*nu;

    w0sq = N*3.1826e3;
    s = 1.7591e11*mag_field;
    w = ang_freq;

    a = (w0sq./nu.^2).*c_integral(w./nu, 3/2);
    b = (5*w0sq./(2*w.*nu)).*c_integral(w./nu, 5/2);
    c = (w0sq.*(w-s)./(w.*nu.^2)).*c_integral((w-s)./nu, 3/2);
    d = (5*w0sq./(2*w.*nu)).*c_integral((w-s)./nu, 5/2);
    e = (w0sq.*(w+s)./(w.*nu.^2)).*c_integral((w+s)./nu, 3/2);
    f = (5*w0sq./(2*w.*nu)).*c_integral((w+s)./nu, 5/2);

    eps1 = 1 - a - 1i*b;
    eps2 = 0.5*(f-d) + 0.5i*(c-e);
    eps3 = a - 0.5*(c+e) + 1i*(b - 0.5*(f+d));

    A = 2*eps1.*(eps1 + eps3);
    B = eps3.*(eps1 + eps3) + eps2.^2;
    C = 2*eps1.*eps2;
    D = 2*eps1;
    E = 2*eps3;

    ref_ind_sq = (A + B.*(sind(theta).^2) + ...
	 OX_mode * sqrt(B.^2 .* sind(theta).^4 - ...
	 (C.*cosd(theta)).^2))./(D + E.*sind(theta).^2);

    chi = abs(imag(sqrt(ref_ind_sq)));

  else    % Appleton-Hartree

    const1 = 80.616385929648814;
    const2 = 2.799249007652821e+10;

    X = const1 * N ./ freq.^2;
    YL = const2 * mag_field .* cosd(theta) ./ freq;
    YT = const2 * mag_field .* sind(theta) ./ freq;
    Z = nu ./ ang_freq;
    iZ = i*Z;
    YT_X_iZ = YT.^2 ./ (2*(1 - X - i*Z));

    ref_ind_sq = 1 - X ./ (1 -iZ - YT_X_iZ + OX_mode*sqrt(YT_X_iZ.^2 + YL.^2));

    chi = abs(imag(sqrt(ref_ind_sq)));

  end

  % The absorption in nepers
  dL = chi.*ds;
  L = sum(dL(~isnan(dL))) * ang_freq / speed_light;

  % Convert absorption to dB
  absorption = L * neper_to_dB;

end

function C = c_integral(x,p)
  de = 0.1;
  e = 0:de:50;

  num = e.^p .* exp(-e);
  f = bsxfun(@rdivide, num, bsxfun(@plus, e.^2, x(:).^2));

  C = ((2/3)*sum(f(:,1:2:end),2) + (4/3)*sum(f(:,2:2:end),2))*de/gamma(p+1);

  C = reshape(C, size(x));
end
