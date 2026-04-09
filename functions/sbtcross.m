function [CrossM] = sbtcross(x,y,es,fs,epsil,psi,Nturns,c,NLegendre,rtolr,atolr)
%SBTCROSS Compute the pairwise SBT interaction matrix between two filaments.
%
%   CrossM = SBTCROSS(x,y,es,fs,epsil,psi,Nturns,c,NLegendre,rtolr,atolr)
%   returns the Legendre-mode slender-body-theory matrix describing the
%   hydrodynamic interaction between two filaments with prescribed
%   positions and orientations.
%
%   INPUTS
%       x           Position of filament 1, size 3-by-1
%       y           Position of filament 2, size 3-by-1
%       es          Orientation matrix of filament 1, size 3-by-3
%       fs          Orientation matrix of filament 2, size 3-by-3
%       epsil       Slenderness parameter
%       psi         Helical angle (set psi = 0 for a straight filament)
%       Nturns      Number of turns (set Nturns = 1 for a straight filament)
%
%   OPTIONAL INPUTS
%       c           Chirality parameter (default: -1)
%       NLegendre   Number of Legendre modes (default: 10)
%       rtolr       Relative tolerance passed to INTEGRAL2 (default: 1e-3)
%       atolr       Absolute tolerance passed to INTEGRAL2 (default: 1e-4)
%
%   OUTPUTS
%       CrossM      Cross-interaction SBT matrix,
%                   size (3*NLegendre)-by-(3*NLegendre)
%
%   NOTES
%       The matrix is assembled by integrating the hydrodynamic kernel
%       against Legendre modes on both filaments.
%
%       The present implementation assumes the two filaments have identical
%       geometric shape parameters.

%%%%%%%%%%%%
% Preamble
%%%%%%%%%%%%
% Check inputs
if nargin<7
    error('Give at least curve shape and filament config.');
end
if nargin<8 || isempty(c)
    c = -1;
end
if nargin<9 || isempty(NLegendre)
    NLegendre = 10;
end
if nargin<10 || isempty(rtolr)
    rtolr = 1e-3;
end
if nargin<11 || isempty(atolr)
    atolr = 1e-4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Legendre polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLegendre = zeros(NLegendre);
MLegendre(1,1) = 1;
MLegendre(2,2) = 1;
for p=3:NLegendre
    MLegendre(p,:) = (2*p-3)/(p-1)*[0 MLegendre(p-1,1:NLegendre-1)] - (p-2)/(p-1)*MLegendre(p-2,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commonly used quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%
piN  = pi*Nturns;
SIN  = sin(psi);
COS  = cos(psi);
R  = SIN/piN;
delta = epsil^2/2;

%%%%%%%%%%%%%%%%%%%%%
% Compute SBT matrix
%%%%%%%%%%%%%%%%%%%%%
% Initialise output
CrossM = zeros(3*NLegendre);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define centreline, tangent, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First filament
r11 = @(s) x(1)+R*cos(piN*s)*es(1,1)+c*R*sin(piN*s)*es(1,2)+COS*s*es(1,3);
r12 = @(s) x(2)+R*cos(piN*s)*es(2,1)+c*R*sin(piN*s)*es(2,2)+COS*s*es(2,3);
r13 = @(s) x(3)+R*cos(piN*s)*es(3,1)+c*R*sin(piN*s)*es(3,2)+COS*s*es(3,3);

% Second filament
r21 = @(s) y(1)+R*cos(piN*s)*fs(1,1)+c*R*sin(piN*s)*fs(1,2)+COS*s*fs(1,3);
r22 = @(s) y(2)+R*cos(piN*s)*fs(2,1)+c*R*sin(piN*s)*fs(2,2)+COS*s*fs(2,3);
r23 = @(s) y(3)+R*cos(piN*s)*fs(3,1)+c*R*sin(piN*s)*fs(3,2)+COS*s*fs(3,3);

% Need norm of R2(s,sd) = r1(s)-r2(sd) (i.e. distance rel to helix 2)
R2n = @(s,sd) sqrt((r11(s)-r21(sd)).^2+(r12(s)-r22(sd)).^2+(r13(s)-r23(sd)).^2);

% Need unit vector \hat{R2} too
R2h1 = @(s,sd) (r11(s)-r21(sd))./R2n(s,sd);
R2h2 = @(s,sd) (r12(s)-r22(sd))./R2n(s,sd);
R2h3 = @(s,sd) (r13(s)-r23(sd))./R2n(s,sd);

tic
for n=0:NLegendre-1
    for m=0:NLegendre-1
        ii = 1; jj = 1;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            ((1+R2h1(s,sd).*R2h1(s,sd))./R2n(s,sd)+...
            delta*(1-3*R2h1(s,sd).*R2h1(s,sd))./R2n(s,sd).^3);
        I = integral2(fun,-1,1,-1,1,'RelTol',rtolr,'AbsTol',atolr);
        
        % Effect on helix 1 due to forces on helix 2
        CrossM(3*m+ii,3*n+jj) = I;
              
        % 22
        ii = 2; jj = 2;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            ((1+R2h2(s,sd).*R2h2(s,sd))./R2n(s,sd)+...
            delta*(1-3*R2h2(s,sd).*R2h2(s,sd))./R2n(s,sd).^3);
        I = integral2(fun,-1,1,-1,1,'RelTol',rtolr,'AbsTol',atolr);
        
        % Effect on helix 1 due to forces on helix 2
        CrossM(3*m+ii,3*n+jj) = I;
        
        % 33
        ii = 3; jj = 3;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            ((1+R2h3(s,sd).*R2h3(s,sd))./R2n(s,sd)+...
            delta*(1-3*R2h3(s,sd).*R2h3(s,sd))./R2n(s,sd).^3);
        I = integral2(fun,-1,1,-1,1,'RelTol',rtolr,'AbsTol',atolr);
        
        % Effect on helix 1 due to forces on helix 2
        CrossM(3*m+ii,3*n+jj) = I;
        
        % 12 and 21
        ii = 1; jj = 2;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            ((0+R2h1(s,sd).*R2h2(s,sd))./R2n(s,sd)+...
            delta*(0-3*R2h1(s,sd).*R2h2(s,sd))./R2n(s,sd).^3);
        I = integral2(fun,-1,1,-1,1,'RelTol',rtolr,'AbsTol',atolr);
        
        % Effect on helix 1 due to forces on helix 2
        CrossM(3*m+ii,3*n+jj) = I;
        CrossM(3*m+jj,3*n+ii) = I;
        
        % 13 and 31
        ii = 1; jj = 3;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            ((0+R2h1(s,sd).*R2h3(s,sd))./R2n(s,sd)+...
            delta*(0-3*R2h1(s,sd).*R2h3(s,sd))./R2n(s,sd).^3);
        I = integral2(fun,-1,1,-1,1,'RelTol',rtolr,'AbsTol',atolr);
        
        % Effect on helix 1 due to forces on helix 2
        CrossM(3*m+ii,3*n+jj) = I;
        CrossM(3*m+jj,3*n+ii) = I;
        
        % 23 and 32
        ii = 2; jj = 3;
        
        % Define the integrand
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            ((0+R2h2(s,sd).*R2h3(s,sd))./R2n(s,sd)+...
            delta*(0-3*R2h2(s,sd).*R2h3(s,sd))./R2n(s,sd).^3);
        I = integral2(fun,-1,1,-1,1,'RelTol',rtolr,'AbsTol',atolr);
        
        % Effect on helix 1 due to forces on helix 2
        CrossM(3*m+ii,3*n+jj) = I;
        CrossM(3*m+jj,3*n+ii) = I;
    end
end
% now = toc;
% disp(['Compute SBT cross-interaction matrix - finished:' num2str(now) 's'])

end

