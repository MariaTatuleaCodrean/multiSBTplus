function [ResM,SbtM,fmom] = sbtself(epsil,psi,Nturns,c,NLegendre,rtolr,atolr)
%SBTSELF Compute single-filament SBT and resistance matrices.
%
%   [ResM,SbtM,fmom,ResMwithmom] = SBTSELF(epsil,psi,Nturns,c,...
%   NLegendre,rtolr,atolr) computes the Legendre-mode slender-body-theory
%   matrix for an isolated filament together with the corresponding
%   rigid-body resistance matrix.
%
%   INPUTS
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
%       ResM        6-by-6 resistance matrix without moment contributions
%       SbtM        SBT matrix in Legendre representation,
%                   size (3*NLegendre)-by-(3*NLegendre)
%       fmom        Force-moment quantity used in asymptotic interaction
%                   calculations
%
%   NOTES
%       The filament centreline is assumed to be either straight or helical.
%
%       The routine first constructs the single-filament SBT matrix and then
%       computes the resistance matrix by solving for the force density
%       associated with each rigid-body translation and rotation mode.

%%%%%%%%%%%%
% Preamble
%%%%%%%%%%%%
% Check inputs
if nargin<3
    error('Give at least curve and slenderness parameter.');
end
if nargin<4 || isempty(c)
    c = -1;
end
if nargin<5 || isempty(NLegendre)
    NLegendre = 10;
end
if nargin<6 || isempty(rtolr)
    rtolr = 1e-3;
end
if nargin<7 || isempty(atolr)
    atolr = 1e-4;
end

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Legendre polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLegendre = zeros(NLegendre);
MLegendre(1,1) = 1;
MLegendre(2,2) = 1;

% Pn = cell(N,1);
% Pn{1} = @(s) s.^0;
% Pn{2} = @(s) s;

for k=3:NLegendre
    MLegendre(k,:) = (2*k-3)/(k-1)*[0 MLegendre(k-1,1:NLegendre-1)] - (k-2)/(k-1)*MLegendre(k-2,:);
    %     Pn{k} = @(s) reshape(MLegendre(k,:)*reshape(s,1,numel(s)).^transpose(0:N-1),size(s));
end

Pnvector = @(s) MLegendre*s.^transpose(0:NLegendre-1);

%%%%%%%%%%%%%%%%
% Define shape
%%%%%%%%%%%%%%%%
piN  = pi*Nturns;
piN2 = piN/2;
SIN  = sin(psi);
SIN2 = SIN^2;
COS  = cos(psi);
COS2 = COS^2;
SINCOS = SIN*COS;
R  = SIN/piN;
R2 = R^2;
RR = 2*R;
RR2 = RR^2;
r1 = @(s) R*cos(piN*s);
r2 = @(s) c*R*sin(piN*s);
r3 = @(s) COS*s;
t1 = @(s) -SIN*sin(piN*s);
t2= @(s) c*SIN*cos(piN*s);
t3 = @(s) COS*s.^0;
titj = @(s) [SIN^2*sin(piN*s)^2,...
    -SIN*sin(piN*s)*c*SIN*cos(piN*s),...
    -SIN*sin(piN*s)*COS;...
    -SIN*sin(piN*s)*c*SIN*cos(piN*s),...
    SIN^2*cos(piN*s)^2,...
    c*SIN*cos(piN*s)*COS;...
    -SIN*sin(piN*s)*COS,...
    c*SIN*cos(piN*s)*COS,...
    COS^2*s^0];

%%%%%%%%%%%%%%%%%%%%%
% Prepare SBT kernel
%%%%%%%%%%%%%%%%%%%%%
% Need norm of R0(s,sd) = r(s)-r(sd)
% R0n = @(s,sd) sqrt((r1(s)-r1(sd)).^2+(r2(s)-r2(sd)).^2+(r3(s)-r3(sd)).^2);
% R0n = @(s,sd) sqrt((R*cos(piN*s)-R*cos(piN*sd)).^2+...
%     (c*R*sin(piN*s)-c*R*sin(piN*sd)).^2+...
%     (COS*s-COS*sd).^2);
% R0n = @(s,sd) sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2);

% Need unit vector \hat{R0} too
% R0h1 = @(s,sd) (R*cos(piN*s)-R*cos(piN*sd))./...
%     sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2);
% R0h2 = @(s,sd) (c*R*sin(piN*s)-c*R*sin(piN*sd))./...
%     sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2);
% R0h3 = @(s,sd) (COS*s-COS*sd)./...
%     sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2);

%%%%%%%%%%%%%%%%%%%%%
% Compute SBT matrix
%%%%%%%%%%%%%%%%%%%%%
% Initialise output
SbtM = zeros(3*NLegendre);

for n=0:NLegendre-1
    % Determine eigenvalue for n
    En = 2*sum(1./(1:n));
    
    % Determine M^3_m,n factor (since it only depends on n)
    M3 = -En+log(4/epsil^2)-3;
    
    % Determine M^2_m,n (since it only depends on n)
    M2 = 2*(M3+4)/(2*n+1);
    
    % Add local contribution
    loc = 3*n+(1:3);
    SbtM(loc,loc) = SbtM(loc,loc) + M2*eye(3);
    
    % Calculate contribution from M^3_m,n,i,j
    fun = @(s) kron(Pnvector(s),titj(s))*(MLegendre(n+1,:)*s.^transpose(0:NLegendre-1));
    SbtM(:,loc) = SbtM(:,loc) + M3*integral(fun,-1,1,'ArrayValued',true);
    
    for m=0:NLegendre-1
        % Keep singularities on the integration boundaries
        %         sdlow = @(s) max(s-0.1*abs(s),-1);
        %         sdupp = @(s) min(s+0.1*abs(s),1);
        sdlim = @(sd) sd;
        
        % 11
        ii = 1; jj = 1;
        
        % Define the integrand
        %         fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:N-1),size(s)).*...
        %             reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:N-1),size(sd)).*...
        %             ((1+R0h1(s,sd).*R0h1(s,sd))./R0n(s,sd)-(1+t1(s).*t1(s))./abs(s-sd));
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (1./sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2)+...
            (R*cos(piN*s)-R*cos(piN*sd)).^2./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (1+SIN2*sin(piN*s).^2)./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        
        % 22
        ii = 2; jj = 2;
        
        % Define the integrand
        %         fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:N-1),size(s)).*...
        %             reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:N-1),size(sd)).*...
        %             ((1+R0h2(s,sd).*R0h2(s,sd))./R0n(s,sd)-(1+t2(s).*t2(s))./abs(s-sd));
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (1./sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2)+...
            (R*sin(piN*s)-R*sin(piN*sd)).^2./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (1+SIN2*cos(piN*s).^2)./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        
        % 33
        ii = 3; jj = 3;
        
        % Define the integrand
        %         fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:N-1),size(s)).*...
        %             reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:N-1),size(sd)).*...
        %             ((1+R0h3(s,sd).*R0h3(s,sd))./R0n(s,sd)-(1+t3(s).*t3(s))./abs(s-sd));
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (1./sqrt(RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2)+...
            (COS*s-COS*sd).^2./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (1+COS2)./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        
        % 12 and 21
        ii = 1; jj = 2;
        
        % Define the integrand
        %         fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:N-1),size(s)).*...
        %             reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:N-1),size(sd)).*...
        %             ((1+R0h1(s,sd).*R0h2(s,sd))./R0n(s,sd)-(1+t1(s).*t2(s))./abs(s-sd));
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (c*R2*(cos(piN*s)-cos(piN*sd)).*(sin(piN*s)-sin(piN*sd))./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (0-c*SIN2*sin(piN*s).*cos(piN*s))./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        I = integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        I = integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        % 13 and 31
        ii = 1; jj = 3;
        
        % Define the integrand
        %         fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:N-1),size(s)).*...
        %             reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:N-1),size(sd)).*...
        %             ((1+R0h1(s,sd).*R0h3(s,sd))./R0n(s,sd)-(1+t1(s).*t3(s))./abs(s-sd));
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (R*COS*(cos(piN*s)-cos(piN*sd)).*(s-sd)./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (0-SINCOS*sin(piN*s))./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        I = integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        I = integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        % 23 and 32
        ii = 2; jj = 3;
        
        % Define the integrand
        %         fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:N-1),size(s)).*...
        %             reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:N-1),size(sd)).*...
        %             ((1+R0h2(s,sd).*R0h3(s,sd))./R0n(s,sd)-(1+t2(s).*t3(s))./abs(s-sd));
        fun = @(s,sd) reshape(MLegendre(m+1,:)*reshape(s,1,[]).^transpose(0:NLegendre-1),size(s)).*...
            reshape(MLegendre(n+1,:)*reshape(sd,1,[]).^transpose(0:NLegendre-1),size(sd)).*...
            (c*R*COS*(sin(piN*s)-sin(piN*sd)).*(s-sd)./...
            (RR2*sin(piN2*(s-sd)).^2+COS2*(s-sd).^2).^(3/2)-...
            (0+c*SINCOS*cos(piN*s))./abs(s-sd));
        
        % Calculate contribution from M^1_m,n,i,j
        I = integral2(fun,-1,1,-1,sdlim,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
        
        I = integral2(fun,-1,1,sdlim,+1,'RelTol',rtolr,'AbsTol',atolr);
        SbtM(3*m+ii,3*n+jj) = SbtM(3*m+ii,3*n+jj) + I;
        SbtM(3*m+jj,3*n+ii) = SbtM(3*m+jj,3*n+ii) + I;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute resistance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise output
ResM = zeros(6);
fmom = zeros(1,6);

for k=1:6
    % Loop over linear velocities and angular velocities in principal
    % directions x,y,z
    U = double(k==[1; 2; 3]);
    Omega = double(k==[4; 5; 6]);
    u = @(s) U + cross(Omega,[r1(s); r2(s); r3(s)],1);
    
    % Calculate velocity components
    fun = @(s) 8*pi*kron(Pnvector(s),u(s));
    a = integral(fun,-1,1,'ArrayValued',true);
    a = reshape(a,3*NLegendre,1);
    
    % Calculate force components
    f = SbtM\a;
    f = reshape(f,3,NLegendre);
    
    % Compute total force
    F = 2*f(:,1);
    
    % Compute total torque
    T = zeros(3,1);
    for n = 1:NLegendre
        fun = @(s) [r1(s); r2(s); r3(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        T = T + cross(integral(fun,-1,1,'ArrayValued',true),f(:,n),1);
    end
       
    if k>=4
        % Add contribution from rotlet singularities
        rot = @(s) 4*pi*epsil^2*(1-s^2)*...
            dot(Omega,[t1(s); t2(s); t3(s)],1)*...
            [t1(s); t2(s); t3(s)];
        T = T + integral(rot,-1,1,'ArrayValued',true);
    end
       
    % Save ResM 
    ResM(1:3,k) = F;
    ResM(4:6,k) = T;

    % Compute force moment F1k
    for n = 1:NLegendre
        fun = @(s) [r1(s); r2(s); r3(s)] .* (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        fmom(k) = fmom(k) + dot(integral(fun,-1,1,'ArrayValued',true),f(:,n),1);
        
        fun = @(s) r1(s) * (MLegendre(n,:)*s.^transpose(0:NLegendre-1));
        fmom(k) = fmom(k) - 3*(integral(fun,-1,1,'ArrayValued',true))*f(1,n);
    end
    
end
now = toc;
disp(['Compute SBT self-resistance matrix - finished:' num2str(now) 's'])
end