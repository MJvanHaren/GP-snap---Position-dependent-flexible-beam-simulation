clear all; close all; clc;
addpath(genpath('functions'));
SetPlotLatexStyle;
%% constants
l = 500; %mm
delta = 10;
n = 6;
n_s = 200;
Ts = 1e-3;
N_trial = 6;
%% inputs
X = linspace(delta,l-delta,n)'; % Training inputs
X_s = linspace(0,l,n_s)'; % Test inputs (for visualization)
%% perform ILCBF multiple time for snap parameter
y = zeros(5,n);
for i = 1:n
        [y(:,i),~] = ILCBFSimscape(X(i),l,Ts,N_trial);
end
%% GP
startup;
Y = y(5,:)';
meanfunc = [];
covfunc = @covSEiso;              % Squared Exponental covariance function
likfunc = @likGauss;              % Gaussian likelihood
hypGuess = struct('mean', [], 'cov', [5e0 -15], 'lik', log(1e0));
hypOpt = minimize(hypGuess, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, X, Y); % optimize hyperparameters
[mu s2] = gp(hypOpt, @infGaussLik, meanfunc, covfunc, likfunc, X, Y, X_s);

figure(2);clf;
f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
fill([X_s; flipdim(X_s,1)], f, [7 7 7]/8)
hold on; plot(X_s, mu,'Linewidth',2); plot(X, Y, '+','Markersize',17,'Linewidth',1.3);

xlabel('Sensor y-position [$mm$]');
ylabel('Snap feedforward parameter [$kg/s^2$]')