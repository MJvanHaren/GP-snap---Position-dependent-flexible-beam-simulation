clear all; close all; clc;
addpath(genpath('functions'));
SetPlotLatexStyle;
%% constants
l = 500; %mm
delta = 10;
n = 6;
n_s = 200;
%% inputs
X = linspace(delta,l-delta,n)'; % Training inputs
X_s = linspace(0,l,n_s)'; % Test inputs (for visualization)
%% perform ILCBF multiple time for snap parameter
y = zeros(5,n);
for i = 1:n
        [y(:,i),~] = ILCBFSimscape(X(i),l);
end
%% GP
