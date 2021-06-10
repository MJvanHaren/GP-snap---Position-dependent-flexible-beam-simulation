clear all; close all; clc;
addpath('functions');
SetPlotLatexStyle;
%% inputs
x0 = 75;
qf = 5e-3;
Ts = 1e-3;

n=30;
N_trial = 3;
beamElements = 1;
[ty,ddy] = make4(qf,1e-3,5e-3,5e-2,2e0,Ts);
[~,t,s,j,a,v,r,~] = profile4(ty,ddy(1),Ts);
% t = 0:Ts:3;
% r=sin(100*t)';
% a = -100^2*sin(100*t)';

% figure
% subplot(2,3,1)
% plot(v);
% subplot(2,3,2)
% plot(a);
% subplot(2,3,3)
% plot(j);
% subplot(2,3,4)
% plot(s);
% subplot(2,3,5)
% plot(r);

N2 = 2^12;% for fft
f = 1/Ts*(0:(N2/2))/N2;
r_fft=fft(r,N2);
P2 = abs(r_fft/N2);
specContent = P2(1:N2/2+1);
specContent(2:end-1) = 2*specContent(2:end-1);
subplot(2,3,6)
semilogx(f,specContent);
grid on; xlabel('Frequency [Hz]');

N = length(t);
Tend = t(end);
%% system
load FBcontroller_09062021
sys = c2d(shapeit_data.P.sys,Ts);
controller = d2d(ss(shapeit_data.C_tf_z),Ts);
PS = feedback(sys,controller,-1);
%% weighting
We          = eye(N)*1e6;
Wf          = eye(N)*1e-10;
WDf         = eye(N)*0e-1;
%% BF

Psi = [r a];
npsi = size(Psi,2);
JPsi = zeros(N,npsi);

for iBF = 1:npsi
    JPsi(:,iBF) = lsim(PS,Psi(:,iBF));
end
R = JPsi.'*We*JPsi+Psi.'*(Wf+WDf)*Psi;
Rinv = eye(size(R,2))/R;

Q = Rinv*(JPsi.'*We*JPsi+Psi.'*WDf*Psi);
L = Rinv*(JPsi.'*We);

%% init ILC with BF
theta_jplus1 = zeros(npsi,1);
f_jplus1 = Psi*theta_jplus1;
%% initialize plotting and storage for ILC
PlotTrialData;

% Initialize storage variables.
history.f           = NaN(N,N_trial);
history.u           = NaN(N,N_trial);
history.e           = NaN(N,N_trial);
history.eNorm       = NaN(1,N_trial);
history.eInfNorm    = NaN(1,N_trial);
history.theta_j     = NaN(npsi,N_trial); 

for trial = 1:N_trial
    f_j = f_jplus1;
    theta_j(:,trial) = theta_jplus1;
    
    out = sim('flexibleBeamILCBF');
    
    % load simulation data:
    u_j = out.simout(:,1);
    
    e_j = out.simout(:,2);
   
    % Store trial data.
    history.f(:,trial)          = f_j;
    history.u(:,trial)          = u_j;
    history.e(:,trial)          = e_j;
    history.eNorm(:,trial)      = norm(e_j,2);
    history.eInfNorm(:,trial)   = norm(e_j,Inf);
    history.theta_j(:,trial)  = theta_j(:,trial);
    
    PlotTrialData;

        theta_jplus1 = (Q*theta_j(:,trial)+L*e_j);
        f_jplus1 = Psi*theta_jplus1;
end


