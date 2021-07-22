clear all;close all; clc;
SetPlotLatexStyle;


L = 0.5;                                                                    % [m], length of free-free beam
W = 40e-3;                                                                  % [m], width (or height) of free-free beam
Th = 1e-3;                                                                  % [m], thickness of free-free beam
Ro = 7850;                                                                  % [Kg/m^3], density of material (steel)                                                                     % [-], shape of beam, 3 means rectangular
mass = L*W*Th*Ro;                                                           % [kg]
s = tf('s');


G0 = ModelFlexibleBeamFirstPrinciple(1);
G125 = ModelFlexibleBeamFirstPrinciple(125);
G250 = ModelFlexibleBeamFirstPrinciple(250);
figure; subplot(211);
coolbode(G0,'mag'); hold on;
coolbode(G125,'mag');coolbode(G250,'mag'); 
% coolbode(1/(mass*s^2),'mag');
legend('$\bar{\rho} = 0$ mm','$\bar{\rho} = 125$ mm','$\bar{\rho} = 250$ mm','mass line')
subplot(212);
coolbode(1/(mass*s^2),'mag'); hold on;
coolbode(G125-1/(mass*s^2),'mag');coolbode(G250-1/(mass*s^2),'mag');
legend('Rigid-Body Mode','Flexible Modes $\bar{\rho}=125$ mm','Flexible Modes $\bar{\rho}=250$ mm')
title('Separation of Rigid-Body and Flexible Modes')
%%
annotation('doublearrow',[0.165 0.165], [0.26 0.37]);
text(1.6,-80,{'Position-Dependent', 'Compliance'},'fontsize',8);