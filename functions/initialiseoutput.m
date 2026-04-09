function [x,es,F,T,U,Om] = initialiseoutput(x0,es0,Nfilaments,Ntsteps)
%INITIALISEOUTPUT Preallocate arrays for filament configuration time histories.
%
%   [x,es,F,T,U,Om] = INITIALISEOUTPUT(x0,es0,Nfilaments,Ntsteps)
%   allocates storage for filament positions, orientations, forces,
%   torques, translational velocities, and angular velocities over a
%   prescribed number of time steps.
%
%   INPUTS
%       x0          Initial filament positions, size 3-by-1-by-Nfilaments
%       es0         Initial filament orientation frames,
%                   size 3-by-3-by-Nfilaments
%       Nfilaments  Number of filaments
%       Ntsteps     Number of stored time steps
%
%   OUTPUTS
%       x           Position history, size 3-by-1-by-Nfilaments-by-Ntsteps
%       es          Orientation history,
%                   size 3-by-3-by-Nfilaments-by-Ntsteps
%       F           Force history, size 3-by-1-by-Nfilaments-by-Ntsteps
%       T           Torque history, size 3-by-1-by-Nfilaments-by-Ntsteps
%       U           Translational velocity history,
%                   size 3-by-1-by-Nfilaments-by-Ntsteps
%       Om          Angular velocity history,
%                   size 3-by-1-by-Nfilaments-by-Ntsteps
%
%   NOTES
%       The initial configuration is stored in the first time slice:
%       x(:,:,:,1) = x0 and es(:,:,:,1) = es0.
%
%       Dimension ordering is:
%           1 - Cartesian component
%           2 - column index within vectors/matrices
%           3 - filament index
%           4 - time index

% Preallocate position and orientation histories.
x = zeros(3,1,Nfilaments,Ntsteps);
es = zeros(3,3,Nfilaments,Ntsteps);

% Store the initial configuration in the first time slice.
x(:,:,:,1) = x0;
es(:,:,:,1) = es0;

% Preallocate histories for force, torque, linear velocity,
% and angular velocity.
F  = zeros(3,1,Nfilaments,Ntsteps); 
T  = zeros(3,1,Nfilaments,Ntsteps);
U  = zeros(3,1,Nfilaments,Ntsteps);
Om = zeros(3,1,Nfilaments,Ntsteps);
end

