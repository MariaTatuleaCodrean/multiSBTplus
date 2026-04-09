function [FullRes] = sbt2res(FullSBT,Nfilaments,alles,psi,Nturns,c,NLegendre)
%SBT2RES Convert a multi-filament SBT matrix into a resistance matrix.
%
%   FullRes = SBT2RES(FullSBT,Nfilaments,alles,psi,Nturns,c,NLegendre)
%   computes the 6N-by-6N resistance matrix corresponding to a full
%   multi-filament slender-body-theory matrix written in Legendre modes.
%
%   INPUTS
%       FullSBT     Full multi-filament SBT matrix,
%                   size (3*NLegendre*Nfilaments)-by-(3*NLegendre*Nfilaments)
%       Nfilaments  Number of filaments
%       alles       Orientation matrices for all filaments,
%                   size 3-by-3-by-Nfilaments
%       psi         Helical angle
%       Nturns      Number of helical turns
%       c           Chirality parameter
%       NLegendre   Number of Legendre modes
%
%   OUTPUTS
%       FullRes     Resistance matrix, size (6*Nfilaments)-by-(6*Nfilaments)
%
%   NOTES
%       The routine computes columns of the resistance matrix by prescribing
%       one rigid-body motion at a time for each filament:
%           - three translational modes
%           - three rotational modes
%
%       Forces and torques are then obtained from the solved Legendre-mode
%       force coefficients.
%
%       This version does not include moment contributions along the
%       filament axis in the final resistance matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Legendre polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLegendre = zeros(NLegendre);
MLegendre(1,1) = 1;
MLegendre(2,2) = 1;
for p=3:NLegendre
    MLegendre(p,:) = (2*p-3)/(p-1)*[0 MLegendre(p-1,1:NLegendre-1)] - (p-2)/(p-1)*MLegendre(p-2,:);
end
Pnvector = @(s) MLegendre*s.^transpose(0:NLegendre-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commonly used quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%
piN  = pi*Nturns;
SIN  = sin(psi);
COS  = cos(psi);
R  = SIN/piN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute resistance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise output
FullRes = zeros(6*Nfilaments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let helix n move, others stationary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Nfilaments
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define centreline, tangent, etc.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    es = alles(:,:,n);

    % First filament (relative to centre of rotation)
    rn1 = @(s) R*cos(piN*s)*es(1,1)+c*R*sin(piN*s)*es(1,2)+COS*s*es(1,3);
    rn2 = @(s) R*cos(piN*s)*es(2,1)+c*R*sin(piN*s)*es(2,2)+COS*s*es(2,3);
    rn3 = @(s) R*cos(piN*s)*es(3,1)+c*R*sin(piN*s)*es(3,2)+COS*s*es(3,3);
       
    for k=1:6
        % Loop over linear velocities and angular velocities in principal
        % directions x,y,z
        U = double(k==[1; 2; 3]);
        Omega = double(k==[4; 5; 6]);
        un = @(s) U + cross(Omega,[rn1(s); rn2(s); rn3(s)],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate velocity components - First, helix n
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fun = @(s) 8*pi*kron(Pnvector(s),un(s));
        an = integral(fun,-1,1,'ArrayValued',true);
        an = reshape(an,3*NLegendre,1);
        
        % Initialise vectors for storing velocity modes
        alla = zeros([size(an) Nfilaments]);
        allm = 1:Nfilaments;
        
        % Save velocity components for helix n
        alla(:,:,n) = an;
              
        % Finished calculation velocity components, put into column form
        a = alla(:);
        
        % Calculate force components
        f = FullSBT\a;
        f = reshape(f,3,NLegendre,Nfilaments); 
        % dim = 1, Cartesian component
        % dim = 2, Legendre mode
        % dim = 3, fil number
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Total force calculation -> resistance matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For every filament
        for m=allm
            % Compute total force
            Fm = 2*f(:,1,m);
            
            % Enter force into resistance matrix
            
            % This is the force exerted by the mth filament, so
            rowloc = (m-1)*6 + (1:3);
            
            % Due to the motion of the nth filament, so
            colloc = (n-1)*6 + k;
            
            % Save in resistance matrix
            FullRes(rowloc,colloc) = Fm;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Total torque calculation -> resistance matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For every filament
        for m=allm
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Define centreline, tangent, etc.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fs = alles(:,:,m);
            
            % Second filament (relative to centre of rotation)
            rm1 = @(s) R*cos(piN*s)*fs(1,1)+c*R*sin(piN*s)*fs(1,2)+COS*s*fs(1,3);
            rm2 = @(s) R*cos(piN*s)*fs(2,1)+c*R*sin(piN*s)*fs(2,2)+COS*s*fs(2,3);
            rm3 = @(s) R*cos(piN*s)*fs(3,1)+c*R*sin(piN*s)*fs(3,2)+COS*s*fs(3,3);
                       
            % Initialise total torque
            Tm = zeros(3,1);
            
            % Compute total torque as sum over Legendre modes
            fm = f(:,:,m);
            for nn = 1:NLegendre
                fun = @(s) [rm1(s); rm2(s); rm3(s)] .* (MLegendre(nn,:)*s.^transpose(0:NLegendre-1));
                Tm = Tm + cross(integral(fun,-1,1,'ArrayValued',true),fm(:,nn),1);
            end
            
            % Note: in v4.0 there was a calculation here for torques due to
            % moments around the filament axis. Removed in v5.0 
            
            % Enter torque into resistance matrix
            
            % This is the torque exerted by the mth filament, so
            rowloc = (m-1)*6 + (4:6);
            
            % Due to the motion of the nth filament, so
            colloc = (n-1)*6 + k;
            
            % Save in resistance matrix
            FullRes(rowloc,colloc) = Tm;
        end
    end
end

end

