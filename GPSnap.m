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
Xtest = [68 248 310 470]';
ntest = length(Xtest);
theta0 = zeros(5,1);
%% perform ILCBF multiple time for snap parameter
y = zeros(5,n);
for i = 1:n
        [y(:,i),e_j] = ILCBFSimscape(X(i),l,Ts,N_trial,theta0);
end
%% GP
startup;

hypGuessMatrix = [5e0 log(mean(abs(y(1,:)))) log(1e-1);
    5e0 log(1e-1) log(1e-2);
    5e0 log(mean(abs(y(3,:)))) log(1e-1);
    5e0 log(mean(abs(y(4,:)))) log(1e-1);
    5e0 -15 log(1e-1)];

for i = 1:5
    Y = y(i,:)';
    if i ==2
        meanfunc = @meanConst;
        meanguess = mean(Y);
    else
        meanfunc = [];
        meanguess = [];
    end
    covfunc = @covSEiso;              % Squared Exponental covariance function
    likfunc = @likGauss;              % Gaussian likelihood
    hypGuess = struct('mean',meanguess, 'cov', hypGuessMatrix(i,1:2), 'lik', hypGuessMatrix(i,3));
    hypOpt(1,i) = minimize(hypGuess, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, X, Y); % optimize hyperparameters
    [mu, s2] = gp(hypOpt(1,i), @infGaussLik, meanfunc, covfunc, likfunc, X, Y, X_s);

    figure(1+i);clf;
    f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
    fill([X_s; flipdim(X_s,1)], f, [7 7 7]/8)
    hold on; plot(X_s, mu,'Linewidth',2); plot(X, Y, '+','Markersize',17,'Linewidth',1.3);

    xlabel('Sensor y-position [$mm$]');
    ylabel('Snap feedforward parameter [$kg/s^2$]')
end
%% estimate on new position
for i = 1:5
    Y = y(i,:)';
    if i ==2
        meanfunc = @meanConst;
        meanguess = mean(Y);
    else
        meanfunc = [];
        meanguess = [];
    end
    [thetaTest(i,:), ~] = gp(hypOpt(1,i), @infGaussLik, meanfunc, covfunc, likfunc, X, Y, Xtest);
%     thetaTest = [repmat(mean(y(1:4,:),2),1,ntest);thetaSnapTest'];
end
%% test on new position
eGP = zeros(length(e_j),ntest);
eConstant = eGP;
eGPSnap = eGP;
eNormGP = zeros(ntest,1);
eNormGPSnap = zeros(ntest,1);
eNormConstant = eNormGP;

for i = 1:ntest
      [~,eGP(:,i)] = ILCBFSimscape(Xtest(i),l,Ts,1,thetaTest(:,i));
      eNormGP(i) = norm(eGP(:,i),2);
      [~,eGPSnap(:,i)] = ILCBFSimscape(Xtest(i),l,Ts,1,[mean(y(1:4,:),2);thetaTest(5,i)]);
      eNormGPSnap(i) = norm(eGPSnap(:,i),2);
      [~,eConstant(:,i)] = ILCBFSimscape(Xtest(i),l,Ts,1,mean(y,2));
      eNormConstant(i) = norm(eConstant(:,i),2);
end
%% visualization
figure
semilogy(eNormGP,'s--','Markersize',15,'Linewidth',1.3)
hold on
semilogy(eNormGPSnap,'o--','Markersize',15,'Linewidth',1.3)
semilogy(eNormConstant,'^--','Markersize',15,'Linewidth',1.3)
xlabel('Test Position Number [-]');
ylabel('$\|e\|_2 \; [m]$');
legend('GP FEEDFORWARD (not only snap!)','GP Snap (not other parameters)','No Position-Dependent feedforward');
