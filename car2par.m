function [a, e, i, OM, om, th] = car2par(rr, vv, mu)

% Trasformation from cartesian coordinates to Keplerian parameters
%
% [a, e, i, OM, om, th] = car2par(rr, vv, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% rr  [ 3x1 ]   position vector            [ km ]
% vv  [ 3x1 ]   velocity vector            [ km/s ]
% mu  [ 1x1 ]   gravitational parameter    [ km^3/s^2 ]
%
%--------------------------------------------------------------------------
% Output arguments:
%
% a   [ 1x1 ]   semi-major axis            [ km ]
% e   [ 1x1 ]   eccentricity               [ - ]
% i   [ 1x1 ]   inclination                [ rad ]
% OM  [ 1x1 ]   RAAN                       [ rad ]
% om  [ 1x1 ]   pericenter anomaly         [ rad ]
% th  [ 1x1 ]   true anomaly               [ rad ]


% Position and speed modules
r = norm(rr);
v = norm(vv);

% Specific mechanical energy and semi-major axis
E = 0.5 * (v^2) - mu / r;
a = (2/r - (v^2) / mu)^(-1);

% Specific angular momentum vector
hh = cross(rr,vv);
h = norm(hh);

% Vector eccentricity and eccentricity
ee = (cross(vv,hh) / mu) - rr/r;
e = norm(ee);

% Inclination
i = acos(hh(3) / h);

% Nodal axis
kk = [0; 0; 1];
NN = (cross(kk,hh) / norm(cross(kk,hh)));

% Longitude of ascending node (RAAN)
if NN(2) >= 0
    OM = acos(NN(1));
else 
    OM = 2*pi - acos(NN(1));
end

% Argument of periapsis
if ee(3) >= 0
    om = acos(dot(NN,ee) / e);
else 
    om = 2*pi - acos(dot(NN,ee) / e);
end

% True anomaly
v_r = dot(vv,rr) / r;

if v_r >= 0
    th = acos(dot(rr,ee) / (r*e));
else
    th = 2*pi - acos(dot(rr,ee) / (r*e));
end

end