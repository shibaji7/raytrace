%
% Name : 
%   ground_bs_loss.m
%
% Purpose :
%   Calculate power loss (dB) of the radio-waves back-scattered from the ground.
%   This is a very simple model of back-scattered loss which is based on two land
%   types viz. land or sea. For the case of sea, the sea is considered to be
%   fully-developed and the back scatter loss is 23dB. For the case of land
%   the back scatter loss is 26dB. An additional loss of 3dB is included to
%   account for polarization mismatch at the receive antenna.
%
%   See : Coleman, "On the simulation of backscatter ionograms", JASTP, vol 59,
%                   2089 - 2099, 1997
%
% Calling sequence :
%   bs_loss = ground_bs_loss(lat, lon);
%
% Inputs :
%   lat      - array of latitudes of start point of rays (deg)  
%   lon      - array of longitudes of start point of rays  (deg)
%
% Outputs :
%   bs_loss - array of the power loss of the radio-waves back-scattered from the 
%             ground (dB)
%
% Author:
%   V1.0  M.A. Cervera  11/09/2006
%
%   V1.1  M.A. Cervera  13/07/2017
%     Allow array inputs
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (ground_bs_loss_matlab_wrapper.for) to the Fortran code (land_type.for).
