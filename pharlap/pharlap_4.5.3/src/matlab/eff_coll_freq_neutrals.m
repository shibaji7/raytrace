%
% Name :
%   eff_coll_freq_neutrals.m
% 
% Purpose / description:
%   Calculates the effective collision frequency between electrons and
%   particular neutral atmospheric species. Uses data from "Effective
%   collision frequency of electrons in gases", Y. Itikawa, Physics of
%   Fluids, vol 16 no. 6, pp 831-835, 1973.
%
% Calling sequence:
%   coll_freq = eff_coll_freq(T_e, number_density, species)
%
% Inputs:
%   T_e             - 1xN array of electron temperatures, in K
%   number_density  - 1xN array of number density of neutral species, in cm^-3
%   species         - a number specifying the neutral species:
%                      1  - N_2, nitrogen
%                      2  - O_2, oxygen
%                      3  - NO, nitric oxide
%                      4  - H_20, water vapour
%                      5  - CO_2, carbon dioxide
%                      6  - CH_4, methane
%                      7  - H, hydrogen
%                      8  - He, helium
%                      9  - O, oxygen
%                      10 - Ar, argon
%
% Output:
%   coll_freq  - the effective electron-neutral collision frequency (Hz)
%                for the input neutral species
%
%   V1.0  L.H. Pederick  31/10/2012
%     Initial version.
%
%   V1.1  M.A. Cervera  18/05/2018
%     Minor code and comment tidy-ups
%

function coll_freq = eff_coll_freq_neutrals(T_e, number_density, species)
    persistent T_e_axis coll_freq_data
    
    if ~isequal(size(T_e), size(number_density))
      error('Inputs 1 and 2 must have the same size')
    end
    
    if (~isnumeric(T_e) || ~isreal(T_e) || ~isnumeric(number_density) || ...
        ~isreal(number_density))
      error('Inputs 1 and 2 must be numeric and real')
    end 
    
    % Electron temperature array for the Itikawa (1973) electron-neutral
    % collision frequency data defined in the coll_freq_data array below.
    if (isempty(T_e_axis))
        T_e_axis = [100  200  300 400 500 1000 1500 2000 2500 3000 3500 ...
	            4000 4500 5000];
    end
    
    % Electron-neutral collision frequency data from Itikawa (1973). Data is 
    % the effective collision frequency as a function of species and electron
    % temerature for unitiary number density of the neutral species.
    if (isempty(coll_freq_data))
        coll_freq_data = [...
            0.255   0.491   0.723   0.948   1.17    2.11    2.86    3.49    4.08    4.72    5.44    6.23    7.08    7.96;   % N_2
            0.123   0.237   0.335   0.425   0.512   0.956   1.41    1.87    2.32    2.72    3.08    3.40    3.68    3.92;   % O_2
            0.589   0.528   0.566   0.673   0.859   2.73    4.88    6.46    7.49    8.13    8.54    8.79    8.95    9.06;   % NO
            158     107     84.9    71.5    62.3    39.3    29.3    23.7    20.1    17.6    15.8    14.5    13.4    12.5;   % H_2O
            10.1    10.1    10.1    9.96    9.76    8.04    6.45    5.34    4.61    4.14    3.86    3.71    3.68    3.74;   % CO_2
            1.51    0.987   0.719   0.568   0.480   0.396   0.486   0.636   0.822   1.04    1.29    1.58    1.91    2.26;   % CH_4         
            3.44    4.85    5.89    6.72    7.42    9.81    11.3    12.2    12.9    13.4    13.8    14.1    14.3    14.5;   % H
            0.448   0.654   0.820   0.963   1.09    1.62    2.05    2.41    2.74    3.03    3.30    3.55    3.78    4.00;   % He
            0.081   0.138   0.192   0.245   0.297   0.551   0.800   1.04    1.28    1.51    1.73    1.94    2.14    2.34;   % O
            0.341   0.293   0.243   0.203   0.173   0.110   0.126   0.185   0.272   0.379   0.504   0.643   0.795   0.960;  % Ar
            ] * 1e-8;
    end
    
    % Calculate the electron-neutral collision frequency for the input
    % neutral species
%    coll_freq = lininterp1f(T_e_axis, coll_freq_data(species, :), T_e, ...
%	                    NaN) .* number_density;
    coll_freq = interp1(T_e_axis, coll_freq_data(species, :), T_e) .* number_density;
end