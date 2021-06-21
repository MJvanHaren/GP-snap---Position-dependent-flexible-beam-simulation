function [G] = ModelFlexibleBeamFirstPrinciple(Ix)
addpath('functions');
% SetPlotLatexStyle;
%% definitions
L = 0.5;                                                                    % [m], length of free-free beam
W = 40e-3;                                                                  % [m], width (or height) of free-free beam
Th = 1e-3;                                                                  % [m], thickness of free-free beam
E = 2.1e11;                                                                 % [Pa], young's modulus of material (steel)
Ro = 7850;                                                                  % [Kg/m^3], density of material (steel)
HMMS = 3;                                                                   % [-], how many mode shapes, amount of modeshapes to estimate using euler equations
CS=3;                                                                       % [-], shape of beam, 3 means rectangular
mass = L*W*Th*Ro;                                                           % [kg]
%% Determine (normalized) modeshapes, damping of modes, freqency of modes and *set* gain matrix P
[Xnx, ~, fn] = FFbeam(L,W,Th,E,Ro,HMMS,CS);
Xnx = Xnx(1:2:end,:);
wn = fn(1:2:end)*2*pi;
Xnxscale = (Xnx'-1)./(max(abs(Xnx'-1)));
factor = sqrt(sum(Xnxscale.^2));

Phi = [-ones(500,1)/sqrt(500) Xnxscale(:,1)./factor(1) Xnxscale(:,2)./factor(2)];
B = [1 zeros(1,498) 1]';
C = [zeros(1,Ix-1) 1 zeros(1,500-Ix)];

M = diag(mass/500*ones(500,1));
Mm = Phi'*M*Phi;
Bm = inv(Mm)*Phi'*B;
Cm = C*Phi;

%% iterating over modes and positions to determine Gy and Gz
s = tf('s');
G = (Cm(1)*Bm(1))/s^2;
n = length(Cm);
zeta = [0.071 0.03];
for r = 2:n
    G = G+2*(Cm(r)*Bm(r))/(s^2+2*zeta(r-1)*wn(r-1)*s+wn(r-1)^2);
end
end