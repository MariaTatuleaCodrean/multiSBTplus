%% demo.m
% DEMO Generate an example input file, run MAIN_CONFIGSTEPPING, and plot results.
%
% This script:
%   1. Generates an input file called 'input_example.mat'
%   2. Runs MAIN_CONFIGSTEPPING on that input
%   3. Loads the resulting resistance matrices
%   4. Plots the first filament configuration using TUBEPLOT
%   5. Plots F_z and T_z on filament 1 versus phase angle phi, assuming
%      all filaments rotate with Omega_z = 1
%
% The geometry used here matches the example values from the original
% input-generation script.
%
% NOTE ON AUTHORSHIP
% Parts of this script were written with assistance from ChatGPT.
% The script was generated from a natural-language description of the 
% desired workflow and from the interfaces of the MATLAB functions provided 
% by the user, which had been developed for the research article cited in 
% the README.md file. Users should review and validate the script before 
% relying on it for scientific use.

clear; clc; close all;
addpath('functions'); 

% Filenames
inputfile = 'input_example.mat';
outputfile = 'output_example.mat';

% Geometric parameters [μm]
r  = 0.2;       % helix radius
p  = 2;         % helix pitch
re = 0.01;      % filament thickness
c  = -1;        % chirality (-1: left-handed, +1: right-handed)
l  = 8;         % filament length
Rb = 0.5;       % bundle radius
M  = 4;         % number of filaments

% Configuration sampling
configresl = 16;

% Optional accuracy parameters
NLegendre = 15;
rtolr = 1e-6;
atolr = 1e-7;

% Generate example input file
makeinputfile(inputfile,r,p,re,c,l,Rb,M,configresl,...
    NLegendre,rtolr,atolr);

%% Run resistance-matrix calculation
main_configstepping(inputfile,outputfile);

%% Load data
load(inputfile,'x','es','psi','Nturns')
load(outputfile,'FullRes')

% Reconstruct sampled phase angles
phi = linspace(0,2*pi,configresl);
phi = phi(1:end-1);
Nconfigs = numel(phi);

%% Plot 1: first configuration in 3D
% Reconstruct the filament centreline from the stored position/orientation
% data and the helical geometry used by the SBT routines.
figure(1); clf;
hold on;

% Non-dimensional helix radius 
R = sin(psi)/(pi*Nturns);

% Discretisation along the filament centreline
s = linspace(-1,1,200);

% Plot each filament in the first stored configuration
for n = 1:M
    x0 = x(:,:,n,1);
    es0 = es(:,:,n,1);

    % Centreline in the lab frame
    ctrline = x0 + ...
        R*cos(pi*Nturns*s).*es0(:,1) + ...
        c*R*sin(pi*Nturns*s).*es0(:,2) + ...
        cos(psi)*s.*es0(:,3);

    % Rescale dimensions by filament half-length
    ctrline = l/2 * ctrline;

    % Plot as a tube with radius set by the filament thickness parameter
    tubeplot(ctrline,re,[0 0.4470 0.7410],24,true);
end

axis equal;
xlabel('x [μm]');
ylabel('y [μm]');
zlabel('z [μm]');               
title('Initial configuration');
view(3);
camlight;
lighting gouraud;
box on;

%% Plot 2: F_z and T_z on filament 1 versus phase angle
% For each configuration, compute the total effect on filament 1 when every
% filament rotates with Omega_z = 1 and all other velocity components are 0.
Fz = zeros(1,Nconfigs+1);
Tz = zeros(1,Nconfigs+1);

for kk = 1:Nconfigs
    Res = FullRes(:,:,kk);

    % Build the generalized velocity vector:
    % each filament contributes [Ux Uy Uz Omegax Omegay Omegaz]
    qdot = zeros(6*M,1);
    for n = 1:M
        qdot(6*(n-1) + 6) = 1;   % Omega_z = 1 for filament n
    end

    % Total forces/torques on all filaments
    FT = Res * qdot;

    % Extract force and torque on filament 1
    Fz(kk) = FT(3);   % z-force on filament 1
    Tz(kk) = FT(6);   % z-torque on filament 1
end

% In SBT code, the lengths were rescaled by the half-filament length.
% When we plot forces and torques, we want to use the full filament length
% L as the chracteristic length scale, to avoid unpleasant factors of 2. 
% Since we prescribe Omega_z = 1 in this demonstration, we have:
Fz = Fz/4;
Tz = Tz/8;

% Append the 2*pi endpoint by periodicity to close the curves.
Fz(end) = Fz(1);
Tz(end) = Tz(1);
phi = [phi phi(1)+2*pi];

figure(2); clf;

subplot(1,2,1)
plot(phi,Fz,'-o','LineWidth',1.5,'MarkerSize',6)
xlabel('\phi')
ylabel('F_z/\mu \Omega L^2')
title('dimensionless axial force')
xlim([0 2*pi])
xticks([0 pi 2*pi])
xticklabels({'0','\pi','2\pi'})
grid on
box on

subplot(1,2,2)
plot(phi,Tz,'-o','LineWidth',1.5,'MarkerSize',6)
xlabel('\phi')
ylabel('T_z/\mu \Omega L^3')
title('dimensionless axial torque')
xlim([0 2*pi])
xticks([0 pi 2*pi])
xticklabels({'0','\pi','2\pi'})
grid on
box on