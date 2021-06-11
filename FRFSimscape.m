clear all; close all; clc;
addpath('functions');
SetPlotLatexStyle;
%% inputs
Ts = 5e-4;
Tend = 5;
t = 0:Ts:Tend;
N = length(t);
l = 500;
POI = 250;
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
%% signals

out = sim('flexibleBeamFRF');
    %%
output = out.simout.Data(:,1:4);
input = out.simout.Data(:,5);
%%
systemIdentification;

%% TF estimate
nfs = 2048*2;
wind = kaiser((N-1)/6,10);
[FRFPOI,~] =tfestimate(input,output(:,4),wind,[],nfs,1/Ts);
[FRFMID,ft] =tfestimate(input,output(:,2),wind,[],nfs,1/Ts);

figure
% semilogx(ft*2*pi,20*log10(abs(FRFPOI))); hold on; grid on;
semilogx(ft*2*pi,20*log10(abs(FRFMID)));
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('POI','Center')
% xlim([1 500])