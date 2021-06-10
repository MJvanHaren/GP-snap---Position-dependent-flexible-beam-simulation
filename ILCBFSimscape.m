function [theta_jplus1,e_j] = ILCBFSimscape(POI,l,Ts,N_trial,theta0)
close all;
%% inputs
if POI < l/2 && POI > 0
    x0 = POI;
elseif POI > l/2 && POI < l
    x0 = POI-l/2;
elseif POI == l/2 || POI == 0 || POI == l
    x0 = l/4;
else
    error('Not good POI measurement!')
end

n=30;

beamElements = 1;
[ty,ddy] = make4(5e-4,1e-3,1e-2,2.5e-1,2e1,Ts); % good choice: 5e-4,1e-3,1e-2,2.5e-1,2e1
[~,t,s,j,a,v,r,~] = profile4(ty,ddy(1),Ts);

% figure
% subplot(2,3,1)
% plot(t,r);
% xlabel('Time [s]');
% ylabel('Reference [$m$]');
% subplot(2,3,2)
% plot(t,v);
% xlabel('Time [s]')
% ylabel('Velocity [$m/s$]')
% subplot(2,3,3)
% plot(t,a);
% xlabel('Time [s]')
% ylabel('Acceleration [$m/s^2$]')
% subplot(2,3,4)
% plot(t,j);
% xlabel('Time [s]')
% ylabel('Jerk [$m/s^3$]')
% subplot(2,3,5)
% plot(t,s);
% xlabel('Time [s]')
% ylabel('Sanp [$m/s^4$]')


N = length(t);
Tend = t(end);

% N2 = 2^16;% for fft
% f = 1/Ts*(0:(N2/2))/N2;
% a_fft=fft(a,N2);
% P2 = abs(a_fft/N2);
% specContent = P2(1:N2/2+1);
% specContent(2:end-1) = 2*specContent(2:end-1);
% subplot(2,3,6)
% semilogx(f,20*log10(specContent));
% grid on; xlabel('Frequency [Hz]');



%% system
load FBcontroller_09062021
sys = c2d(shapeit_data.P.sys,Ts);
controller = d2d(ss(shapeit_data.C_tf_z),Ts);
PS = feedback(sys,controller,-1);
%% weighting
We          = eye(N)*1e6;
Wf          = eye(N)*0e-10;
WDf         = eye(N)*0e-1;
%% BF

Psi = [r v a j s];
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
theta_jplus1 = theta0;
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
%     if x0<250
%         out = sim('flexibleBeamILCBFPOIL','SrcWorkspace','current');
%     elseif x0> 250 && x0 < 500
%         out = sim('flexibleBeamILCBFPOIR','SrcWorkspace','current');
%     else
%         error('These beam lenghts are incorrect or WIP (0/250/500)')
%     end
        out = sim('flexibleBeamILCBF','SrcWorkspace','current');
    
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
end


