function [rotatedSbtM] = rotateSBT(SbtM,es)
%ROTATESBT Rotate a filament SBT matrix from body coordinates to lab coordinates.
%
%   rotatedSbtM = ROTATESBT(SbtM,es) transforms the single-filament
%   slender-body-theory matrix SbtM from a body-fixed frame into the lab
%   frame defined by the orientation matrix es.
%
%   INPUTS
%       SbtM        SBT matrix in the body-fixed frame,
%                   size (3N)-by-(3N)
%       es          Rotation/orientation matrix, size 3-by-3
%
%   OUTPUTS
%       rotatedSbtM Rotated SBT matrix in the lab frame,
%                   size (3N)-by-(3N)
%
%   NOTES
%       The transformation is applied blockwise using Kronecker products.
%       A warning is issued if es is not approximately orthogonal with
%       determinant +1.

% Check that es is approximately a proper rotation matrix.
if abs(det(es)-1) > 1e-3
   warning('Body-frame matrix is no longer orthogonal.') 
end

% Number of Legendre/vector blocks in the SBT matrix.
N = round(size(SbtM,1)/3);

% Build block rotation operators acting on all 3-by-1 vector blocks.
backrotate = kron(eye(N),es);
rotate = kron(eye(N),inv(es));

% Apply the change of basis from body frame to lab frame.
rotatedSbtM = backrotate*SbtM*rotate;

end

