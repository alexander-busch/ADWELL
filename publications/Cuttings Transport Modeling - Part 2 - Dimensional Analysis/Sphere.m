clear all;
clc;
addpath('E:\OneDrive_NTNU\SWwd\MATLAB\Generic\m-files\dimension');

% Define the variable names and the respective dimensions
N = {'D', 'v', 'nu', 'rho', 'd'};
u = {'N', 'm/s', 'm2/s', 'kg/m3', 'm'};

N = {'D', 'v', 'eta', 'rho', 'd'};
u = {'N', 'm/s', 'Pas', 'kg/m3', 'm'};

% create dimensional matrix and conversion factors vector
% order of the dimensions
%   1: mass
%   2: length
%   3: time
%   4: temperatur
%   5: current
%   6: light power
%   7: amount of substance
[d,f] = unit2si(u);

% create the relevance list
% (structure with varable names, base dimensions, conversion factors)
RL = rlist(N,d,f);

% choose the base variables
bv = {'v','d','rho'};

% do the dimensional analysis
piset = diman(RL,bv);

% pretty print the pis
pretty(piset);