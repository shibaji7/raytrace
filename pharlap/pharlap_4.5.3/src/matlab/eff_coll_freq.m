%
% Name :
%   eff_coll_freq.m
% 
% Purpose / description:
%   Calculates the effective collision frequency between electrons and
%   atmospheric constituents (both ions and neutrals).
%
% Calling sequence:
%   coll_freq = eff_coll_freq(T_e, T_ion, elec_dens, neutral_dens)
%
% Inputs:
%   T_e           - 1xN array of electron temperatures (K)
%   T_ion         - 1xN array of ion temperatures (K)
%   neutral_dens  - MxN array of number density of neutral species, in cm^-3
%                     M = 1 : Helium (He)
%                     M = 2 : Atomic Oxygen (O)
%                     M = 3 : Nitrogen (N2)
%                     M = 4 : Oxygen (O2)
%                     M = 5 : Argon
%                     M = 6 : Unused
%                     M = 7 : Atomic Hydrogen (H)
%                   The routine nrlmsise00 will generate this array.
%
% Output:
%   coll_freq  - the effective electron-neutral collision frequency (Hz)
%                for the input neutral species
%
%   V1.0  M.A. Cervera  10/09/2018
%     Initial version
%
%   V1.1  M.A. Cervera  30/07/2019
%     Minor update to allow routine to return NaN for collision frequency
%     if the input electron density is zero
%

function coll_freq = eff_coll_freq(T_e, T_ion, elec_dens, neutral_dens)

  % effective collision frequency between electrons and various
  nu_eN2 = eff_coll_freq_neutrals(T_e, neutral_dens(3,:), 1);
  nu_eO2 = eff_coll_freq_neutrals(T_e, neutral_dens(4,:), 2);
  nu_eH  = eff_coll_freq_neutrals(T_e, neutral_dens(7,:), 7);
  nu_eHe = eff_coll_freq_neutrals(T_e, neutral_dens(1,:), 8);
  nu_eO  = eff_coll_freq_neutrals(T_e, neutral_dens(2,:), 9);
  nu_eAr = eff_coll_freq_neutrals(T_e, neutral_dens(5,:), 10);
  
  % Effective electron-ion collision frequency. 
  nu_ei = eff_coll_freq_ion(T_e, T_ion, elec_dens);

  % total effective electron collision frequency
  coll_freq = nu_eN2 + nu_eO2 + nu_eO + nu_eH + nu_eHe + nu_eAr + nu_ei;      

end