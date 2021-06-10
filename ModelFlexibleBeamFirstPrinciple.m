function [G] = ModelFlexibleBeamFirstPrinciple(Ix)
addpath('functions');
SetPlotLatexStyle;
%% definitions
L = 0.5;                                                                    % [m], length of free-free beam
W = 40e-3;                                                                  % [m], width (or height) of free-free beam
Th = 2e-3;                                                                  % [m], thickness of free-free beam
E = 2.1e11;                                                                 % [Pa], young's modulus of material (steel)
Ro = 7850;                                                                  % [Kg/m^3], density of material (steel)
HMMS = 3;                                                                   % [-], how many mode shapes, amount of modeshapes to estimate using euler equations
CS=3;                                                                       % [-], shape of beam, 3 means rectangular
mass = L*W*Th*Ro;                                                           % [kg]
%% Determine (normalized) modeshapes, damping of modes, freqency of modes and *set* gain matrix P
[Xnx, betaN, fn] = FFbeam(L,W,Th,E,Ro,HMMS,CS);
%% iterating over modes and positions to determine Gy and Gz
s = tf('s');
% figure
% Ix = [125 250 375 500];                                                           % indices on beam for Gy and Gz estimation. 1 means beginning of beam, 721 is the end
% for i = 1:length(Ix)
%     G{i} = 1/(mass*s^2);
%     for r = 1:2:HMMS
%         G{i} = G{i}+((Xnx(r,Ix(i))-1))/(mass*(s^2+2*betaN(r)*s)+(fn(r)*2*pi)^2);
%     end
%     bode(G{i});hold on;
% end
G = 1/(mass*s^2);
for r = 1:2:HMMS
    G = G+((Xnx(r,Ix)-1))/(mass*(s^2+2*betaN(r)*s)+(fn(r)*2*pi)^2);
end
end