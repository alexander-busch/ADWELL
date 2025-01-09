clear all;
clc;
addpath('E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files\dimension');

%% Sphere drag example
VarNam = {'D', 'v', 'nu', 'rho', 'd_p'};
VarDim = {'N', 'm/s', 'm2/s', 'kg/m3', 'm'};

VarNam = {'D', 'v', 'mu', 'rho', 'd_p'};
VarDim = {'N', 'm/s', 'Pas', 'kg/m3', 'm'};

%% Ozbayoglu et al. (2007) 
VarNam = {'A_ratio', 'Q', 'theta', 'D_o', 'D_i', 'rho_f', 'eta_f', 'rpm', 'rho_p', 'd_p', 'g'};
VarDim = {'1', 'm3/s', 'rad', 'm', 'm', 'kg/m3', 'Pas', '1/s', 'kg/m3', 'm', 'GU'};

VarNam = {'A_ratio', 'U', 'theta', 'D_o', 'D_i', 'rho_f', 'eta_f', 'rpm', 'rho_p', 'd_p', 'g'};
VarDim = {'1', 'm/s', 'rad', 'm', 'm', 'kg/m3', 'Pas', '1/s', 'kg/m3', 'm', 'GU'};

%% Cuttings transport
VarNam = {'L', 'e_y', 'e_z', 'm_dot_f', 'rho_f', 'mu_f', 'm_dot_p', 'rho_p', 'd_p', 'd_i', 'd_o', 'phi', 'rpm', 'ROP', 'g'};
VarDim = {'m', '1', '1', 'kg/s', 'kg/m3', 'Pas', 'kg/s', 'kg/m3', 'm', 'm', 'm', 'rad', '1/s', 'm/s', 'GU'};

% Herschel bulkley
VarNam = {'d_h', 'e_y', 'e_z', 'phi', 'U_f', 'rho_f', 'eta_0', 'gamma_dot_0', 'Bi', 'n', 'd_s', 'rho_p-rho_f', 'U_r', 'U_s', 'g', 'alpha_b'};
VarDim = {'m', 'm', 'm', 'rad', 'm/s', 'kg/m3', 'Pas', '1/s', '1', '1', 'm', 'kg/m3', 'm/s', 'm/s', 'GU', '1'};
% Define base variables
bv = {'rho_f','U_f','d_h'};


% create dimensional matrix and conversion factors vector
% order of the dimensions
%   1: mass
%   2: length
%   3: time
%   4: temperatur
%   5: current
%   6: light power
%   7: amount of substance
[d,f] = unit2si(VarDim);

% create the relevance list
% (structure with varable names, base dimensions, conversion factors)
RL = rlist(VarNam,d,f);

% do the dimensional analysis
piset = diman(RL,bv);

%Switch off warnings
%[a, MSGID] = lastwarn();
MSGID = 'symbolic:sym:sym:DeprecateExpressions';
warning('off', MSGID);

% pretty print the pis
pretty(piset);

A1=piset.B;
A2=piset.A;
B=piset.D;
C=piset.C;
C_trans = C';

 p1{ii} = [p1{ii} '*' piset.Name{jj} '^(' num2str(DC(ii,jj)) ')'];