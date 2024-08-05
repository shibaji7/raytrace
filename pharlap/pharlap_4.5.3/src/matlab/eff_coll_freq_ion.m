%
% Name :
%   eff_coll_freq_ion.m
% 
% Purpose / description:
%   Calculates the effective collision frequency between electrons and ions in 
%   the Earth's Ionosphere. Based on Schunk, R. W., and A. F. Nagy (1978), 
%   Electron temperatures in the F region of the ionosphere: Theory and 
%   observations, Rev. Geophys. Space. Phys., 16(3), 355–399.
%
% Calling sequence:
%   coll_freq_ion = eff_coll_freq_ion(T_e, T_ion, elec_dens)
%
% Inputs:
%   T_e        - array of electron temperatures (K)
%   T_ion      - array of ion temperatures (K)
%   elec_dens  - array of number density of electrons (m^-3)
%
% Output:
%   coll_freq  - The effective electron-ion collision frequency (Hz)
%
%   V1.0  M.A. Cervera  18/05/2018
%     Initial version.
%

function coll_freq_ion = eff_coll_freq_ion(T_e, T_ion, elec_dens)

  if ~isequal(size(T_e), size(elec_dens))
    error('All inputs must have the same size')
  end
    
  if ~isequal(size(T_e), size(T_ion))
    error('All inputs must have the same size')
  end

  if (~isnumeric(T_e) || ~isreal(T_e) || ~isnumeric(T_ion) || ...
      ~isreal(T_ion) || ~isnumeric(elec_dens) || ~isreal(elec_dens))
    error('All inputs must be numeric and real')
  end 

  ki_sq = 2.09985255e-4 * elec_dens ./ T_ion;
  ke_sq = 2.09985255e-4 * elec_dens ./ T_e;
  
  ln_coulomb_int = 13.484870477617616 + log(T_e) - 0.5*log(ke_sq) - ...
          ((ke_sq + ki_sq) ./ ki_sq) .* (0.5 * log((ki_sq + ke_sq) ./ ke_sq));
	       
  coll_freq_ion = 3.63315e-6 * elec_dens .* (T_e.^(-3/2)) .* ln_coulomb_int;

end
