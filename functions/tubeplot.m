function [X,Y,Z] = tubeplot(ctrline,r,col,n,ploton)
%TUBEPLOT Generate and optionally plot a tubular surface around a centreline.
%
%   [X,Y,Z] = TUBEPLOT(ctrline,r,col,n,ploton) constructs a surface mesh
%   for a tube of radius r around a discretised 3-D centreline. If PLOTON
%   is true, the surface is also drawn using SURF.
%
%   INPUTS
%       ctrline     Discretised centreline, size 3-by-N
%       r           Tube radius
%
%   OPTIONAL INPUTS
%       col         RGB colour vector for plotting (default: [1 0 0])
%       n           Number of points around each cross-section (default: 16)
%       ploton      Logical flag controlling plotting (default: true)
%
%   OUTPUTS
%       X           Surface x-coordinates
%       Y           Surface y-coordinates
%       Z           Surface z-coordinates
%
%   NOTES
%       The centreline is first resampled using spline interpolation, then
%       local tangent, normal, and binormal directions are used to build
%       the tube surface.
%
%       Hemispherical end caps are added at both ends of the tube.
%
%   EXAMPLE
%       [X,Y,Z] = tubeplot(ctrline,0.05,[0 0 1],24,true);

if nargin < 2 || isempty(ctrline) || isempty(r)
    error('Must give at least centreline and filament thickness.')
else
    if nargin < 5 || isempty(ploton)
        ploton = true; % assume we do plot, unless told otherwise
        if nargin < 4 || isempty(n)
            n = 16;
            if nargin < 3 || isempty(col)
                col = [1 0 0];
            end
        end
    end
end

% Break centreline into coordinates
x = ctrline(1,:);
y = ctrline(2,:);
z = ctrline(3,:);

% Find arclength
sdiff = sqrt((x(2:end)-x(1:end-1)).^2+...
    (y(2:end)-y(1:end-1)).^2+...
    (z(2:end)-z(1:end-1)).^2);
sdiff = [0, sdiff];
s = cumsum(sdiff);

% Interpolate and sample points evenly
N = min(length(x),1000);
ss = linspace(s(1),s(end),N);
xx = spline(s,x,ss);
yy = spline(s,y,ss);
zz = spline(s,z,ss);

% Find tangent
t = xx(3:end)-xx(1:end-2);
u = yy(3:end)-yy(1:end-2);
v = zz(3:end)-zz(1:end-2);

t = [xx(2)-xx(1), t, xx(end)-xx(end-1)];
u = [yy(2)-yy(1), u, yy(end)-yy(end-1)];
v = [zz(2)-zz(1), v, zz(end)-zz(end-1)];

normt = sqrt(t.^2+u.^2+v.^2);
t = t./normt;
u = u./normt;
v = v./normt;

tangent = [t;u;v];

% Find normal
m = t(3:end)-t(1:end-2);
p = u(3:end)-u(1:end-2);
q = v(3:end)-v(1:end-2);

normn = sqrt(m.^2+p.^2+q.^2);
m = m./normn;
p = p./normn;
q = q./normn;

normal = [m;p;q];

% Find binormal
binormal = cross(tangent(:,2:end-1),normal);
b = binormal(1,:);
c = binormal(2,:);
d = binormal(3,:);

% Calculate points in mid-section
alpha = 2*pi*(0:n-1)'/(n-1);
X = xx(2:end-1)+r*(cos(alpha)*m+sin(alpha)*b);
Y = yy(2:end-1)+r*(cos(alpha)*p+sin(alpha)*c);
Z = zz(2:end-1)+r*(cos(alpha)*q+sin(alpha)*d);

% Add smooth caps

% Start of filament
beta = pi/2*(0:4)/4;
% A) flat end that goes exactly through given endpoint
% r0 = sqrt((xx(2)-xx(1))^2+(yy(2)-yy(1))^2+(zz(2)-zz(1))^2);
% B) rounded (hemispherical) cap
r0 = r;
X0 = xx(2)-t(1)*r0*cos(beta)+r*(cos(alpha)*m(1)+sin(alpha)*b(1))*sin(beta);
Y0 = yy(2)-u(1)*r0*cos(beta)+r*(cos(alpha)*p(1)+sin(alpha)*c(1))*sin(beta);
Z0 = zz(2)-v(1)*r0*cos(beta)+r*(cos(alpha)*q(1)+sin(alpha)*d(1))*sin(beta);

% End of filament
beta = pi/2*(4:-1:0)/4;
% A) flat end that goes exactly through given endpoint
% rend = sqrt((xx(end)-xx(end-1))^2+(yy(end)-yy(end-1))^2+(zz(end)-zz(end-1))^2);

% B) rounded cap (hemispherical) 
rend = r;
Xend = xx(end-1)+t(end)*rend*cos(beta)+r*(cos(alpha)*m(end)+sin(alpha)*b(end))*sin(beta);
Yend = yy(end-1)+u(end)*rend*cos(beta)+r*(cos(alpha)*p(end)+sin(alpha)*c(end))*sin(beta);
Zend = zz(end-1)+v(end)*rend*cos(beta)+r*(cos(alpha)*q(end)+sin(alpha)*d(end))*sin(beta);

% Put together mesh
X = [X0 X Xend]; 
Y = [Y0 Y Yend];
Z = [Z0 Z Zend];

% Plot surface
C(:,:,1) = ones(size(X))*col(1);
C(:,:,2) = ones(size(X))*col(2);
C(:,:,3) = ones(size(X))*col(3);

if ploton
    surf(X,Y,Z,C,'EdgeColor','none')
end
end

