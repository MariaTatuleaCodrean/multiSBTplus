function makeinputfile(inputfilename,r,p,re,c,l,Rb,M,configresl,NLegendre,rtolr,atolr)
%MAKEINPUTFILE Generate an input MAT-file for MAIN_CONFIGSTEPPING.
%
%   MAKEINPUTFILE(inputfilename,r,p,re,c,l,Rb,M,configresl)
%   generates a MAT-file containing the geometric, numerical, and
%   configuration data required by MAIN_CONFIGSTEPPING.
%
%   MAKEINPUTFILE(...,NLegendre,rtolr,atolr) additionally
%   specifies the numerical accuracy parameters used when precomputing the
%   isolated-filament SBT quantities.
%
%   INPUTS
%       inputfilename   Name of the MAT-file to create
%       r               Helix radius
%       p               Helix pitch
%       re              Filament thickness
%       c               Chirality
%       l               Filament length
%       Rb              Bundle radius
%       M               Number of filaments
%       configresl      Number of sampled phase values over [0,2*pi]
%
%   OPTIONAL INPUTS
%       NLegendre       Number of Legendre modes (default: 15)
%       rtolr           Relative tolerance for INTEGRAL2 (default: 1e-6)
%       atolr           Absolute tolerance for INTEGRAL2 (default: 1e-7)
%
%   OUTPUTS
%       No direct output arguments.
%       Data are written to INPUTFILENAME for later use by
%       MAIN_CONFIGSTEPPING.
%
%   SAVED VARIABLES
%       status, rtolr, atolr, NLegendre, epsil, psi, Nturns, c,
%       ResM, SbtM, fmom, x, es
%
%   NOTES
%       This function places M identical helical filaments uniformly on a
%       circle of radius Rb and generates a sequence of configurations in
%       which all filaments share the same phase angle.
%
%       Positions are saved with size 3-by-1-by-M-by-Nconfigs and
%       orientations with size 3-by-3-by-M-by-Nconfigs, consistent with
%       MAIN_CONFIGSTEPPING and the other functions in this codebase.
%
%   EXAMPLE
%       make_configstepping_input('input.mat',0.39/2,2.22,0.012,-1,8,0.4,6,16)
%
%       make_configstepping_input('input.mat',0.39/2,2.22,0.012,-1,8,0.4,6,16,...
%           15,1e-6,1e-7)

% Validate required inputs.
if nargin < 9
    error(['Not enough input arguments. Required inputs are: ' ...
        'inputfilename, r, p, re, c, l, Rb, M, configresl.']);
end

% Set optional accuracy parameters to script defaults if not provided.
if nargin < 10 || isempty(NLegendre)
    NLegendre = 15;
end
if nargin < 11 || isempty(rtolr)
    rtolr = 1e-6;
end
if nargin < 12 || isempty(atolr)
    atolr = 1e-7;
end

% Mark the input file as not yet processed.
status = 'notstarted';

% -------------------------------------------------------------------------
% Convert geometric parameters to SBT inputs
% -------------------------------------------------------------------------

Nturns = l/sqrt((2*pi*r)^2 + p^2);
epsil  = 2*re/l;
psi    = atan(2*pi*r/p);
RB     = 2*Rb/l;

% Note: lengths are non-dimensionalised by filament half-length (l/2)

% -------------------------------------------------------------------------
% Precompute isolated-filament SBT quantities
% -------------------------------------------------------------------------

[ResM,SbtM,fmom] = sbtself(epsil,psi,Nturns,c,NLegendre,rtolr,atolr);

% -------------------------------------------------------------------------
% Define the sequence of configurations
% -------------------------------------------------------------------------

% Sample phase angles over [0,2*pi), excluding the repeated endpoint.
phi = linspace(0,2*pi,configresl);
phi = phi(1:end-1);

% Number of stored configurations.
Nconfigs = length(phi);

% -------------------------------------------------------------------------
% Filament positions
% -------------------------------------------------------------------------

% Place filament centres uniformly on a circle in the xy-plane.
x0 = zeros(3,1,M);
for n = 1:M
    theta = 2*pi*n/M;
    x0(:,1,n) = RB * [cos(theta); sin(theta); 0];
end

% Replicate the same positions for all configurations.
x = repmat(x0,1,1,1,Nconfigs);

% -------------------------------------------------------------------------
% Filament orientations
% -------------------------------------------------------------------------

% Store one body-frame orientation matrix for each filament and
% configuration. Here all filaments are taken to be in phase.
es = zeros(3,3,M,Nconfigs);
for mm = 1:Nconfigs
    phi1 = phi(mm);

    % Rotation about the lab-frame z-axis.
    Rz = [cos(phi1), -sin(phi1), 0; ...
          sin(phi1),  cos(phi1), 0; ...
                 0 ,         0 , 1];

    for n = 1:M
        es(:,:,n,mm) = Rz;
    end
end

% -------------------------------------------------------------------------
% Save all variables needed by MAIN_CONFIGSTEPPING
% -------------------------------------------------------------------------

save(inputfilename,'status','epsil','psi','Nturns','c','NLegendre', ...
    'rtolr','atolr','x','es','ResM','SbtM','fmom')
end