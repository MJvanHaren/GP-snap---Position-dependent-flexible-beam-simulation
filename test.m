clear all
close all
G0 = ModelFlexibleBeamFirstPrinciple(1);
G125 = ModelFlexibleBeamFirstPrinciple(125);
G250 = ModelFlexibleBeamFirstPrinciple(250);
figure; subplot(211);
coolbode(G0,'mag'); hold on;
coolbode(G125,'mag');coolbode(G250,'mag');
subplot(212);
s = tf('s');
L = 0.5;                                                                    % [m], length of free-free beam
W = 40e-3;                                                                  % [m], width (or height) of free-free beam
Th = 1e-3;                                                                  % [m], thickness of free-free beam
E = 2.1e11;                                                                 % [Pa], young's modulus of material (steel)
Ro = 7850;                                                                  % [Kg/m^3], density of material (steel)
HMMS = 3;                                                                   % [-], how many mode shapes, amount of modeshapes to estimate using euler equations
CS=3;                                                                       % [-], shape of beam, 3 means rectangular
mass = L*W*Th*Ro;                                                           % [kg]
coolbode(112.8349/(s^2),'mag'); hold on;
coolbode(G125-112.8349/(s^2),'mag');coolbode(G250-112.8349/(s^2),'mag');