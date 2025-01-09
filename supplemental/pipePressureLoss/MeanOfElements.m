function [ out ] = MeanOfElements( v, rep )
%MeanOfElements Summary of this function goes here
%   Detailed explanation goes here

n = numel(v); % Number of elements in x
m = n/rep; % Number of elements in subset of x to average
out = mean(vec2mat(v,n/rep))'; % 

end
