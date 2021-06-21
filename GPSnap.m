clear all; close all; clc;
addpath(genpath('functions'));
SetPlotLatexStyle;
%% constants
l = 500; %mm
delta = 10; % plz >1
n = 5; % preferably (must be) an uneven number so center will be identified
n_s = 200;
Ts = 1e-3;
N_trial = 6;
%% inputs
X = round(linspace(delta,l-delta,n)'); % Training inputs
X_s = linspace(0,l,n_s)'; % Test inputs (for visualization)
Xtest = [30 110 180 248 380 445]'; % mm
ntest = length(Xtest);
%% ref + basis
[ty,ddy] = make4(5e-4,1e-3,1e-2,2.5e-1,2e1,Ts); % good choice: 5e-4,1e-3,1e-2,2.5e-1,2e1
[~,t,s,j,a,v,r,~] = profile4(ty,ddy(1),Ts);
Psi = [a s];
npsi = size(Psi,2);
theta0 = zeros(npsi,1);
if true
    figure
    subplot(2,3,1)
    plot(t,r);
    xlabel('Time [s]');ylabel('Reference [$m$]');
    subplot(2,3,2)
    plot(t,v);
    xlabel('Time [s]');ylabel('Velocity [$m/s$]');
    subplot(2,3,3);
    plot(t,a);
    xlabel('Time [s]');ylabel('Acceleration [$m/s^2$]');
    subplot(2,3,4)
    plot(t,j);
    xlabel('Time [s]');ylabel('Jerk [$m/s^3$]');
    subplot(2,3,5)
    plot(t,s);
    xlabel('Time [s]');ylabel('Sanp [$m/s^4$]');
    N2 = 2^16;% for fft
    f = 1/Ts*(0:(N2/2))/N2;
    a_fft=fft(a,N2);
    P2 = abs(a_fft/N2);
    specContent = P2(1:N2/2+1);
    specContent(2:end-1) = 2*specContent(2:end-1);
    subplot(2,3,6)
    semilogx(f,20*log10(specContent));
    grid on; xlabel('Frequency [Hz]');
end
%% perform ILCBF multiple time for snap parameter
y = zeros(npsi,n);
for i = 1:n
    [y(:,i),e_j] = ILCBFSimscape(X(i),l,Ts,N_trial,theta0,r,Psi,t);
end
%% GP
startup;

Y = y(end,:)';
meanfunc = [];
covfunc = @covSEiso;              % Squared Exponental covariance function
likfunc = @likGauss;              % Gaussian likelihood
hypGuess = struct('mean',[], 'cov', [5e0 log(mean(abs(Y)))], 'lik', log(1e-9));
hypOpt = minimize(hypGuess, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, X, Y); % optimize hyperparameters
[mu, s2] = gp(hypOpt, @infGaussLik, meanfunc, covfunc, likfunc, X, Y, X_s);

figure(2);clf;
f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
fill([X_s; flipdim(X_s,1)], f, [7 7 7]/8)
hold on; plot(X_s, mu,'Linewidth',2); plot(X, Y, '+','Markersize',17,'Linewidth',1.3);

xlabel('Sensor y-position [$mm$]');
ylabel('Snap feedforward parameter [$kg/s^2$]')

%% estimate on new position

[thetaSnapTest, ~] = gp(hypOpt, @infGaussLik, meanfunc, covfunc, likfunc, X, Y, Xtest);
thetaTest = [repmat(y(1:end-1,find(X==l/2)),1,ntest);thetaSnapTest'];
%% test on new position
eGP = zeros(length(e_j),ntest);
eConstant = eGP;
eNormGP = zeros(ntest,1);
eNormConstant = eNormGP;

for i = 1:ntest
    [~,eGP(:,i)] = ILCBFSimscape(Xtest(i),l,Ts,1,thetaTest(:,i),r,Psi,t);
    eNormGP(i) = norm(eGP(:,i),2);
    [~,eConstant(:,i)] = ILCBFSimscape(Xtest(i),l,Ts,1,y(:,find(X==l/2)),r,Psi,t);
    eNormConstant(i) = norm(eConstant(:,i),2);
end
%% visualization
figure(3);clf;
semilogy(Xtest,eNormGP,'s--','Markersize',15,'Linewidth',1.3)
hold on
semilogy(Xtest,eNormConstant,'^--','Markersize',15,'Linewidth',1.3)
xlabel('Scheduling Variable $\rho$ [$mm$]');
ylabel('$\|e\|_2$ [$m$]');
legend('GP Snap Feefdorward','Position-Independent Feedforward');
% ylim([3.8e-6 8.92e-6])
