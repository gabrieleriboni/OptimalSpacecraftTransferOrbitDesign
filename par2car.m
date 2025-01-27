function [rr, vv] = par2car(a, e, i, OM, om, th, mu)

% Trasformation from  Keplerian parameters to cartesian coordinates
%
% [rr, vv] = car2par(a, e, i, OM, om, th)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a   [ 1x1 ]   semi-major axis           [ km ]
% e   [ 1x1 ]   eccentricity              [ - ]
% i   [ 1x1 ]   inclination               [ rad ]
% OM  [ 1x1 ]   RAAN                      [ rad ]
% om  [ 1x1 ]   pericenter anomaly        [ rad ]
% th  [ 1x1 ]   true anomaly              [ rad ]
%
%--------------------------------------------------------------------------
% Output arguments:
%
% rr  [ 3x1 ]   position vector            [ km ]
% vv  [ 3x1 ]   velocity vector            [ km/s ]
% mu  [ 1x1 ]   gravitational parameter    [ km^3/s^2 ]


% Semi-latus rectum
p = a * (1 - (e^2));

% Position module
r = p / (1 + e*cos(th));

% Vectors r_pf, v_pf
rr_pf = r * [cos(th); sin(th); 0];
vv_pf = sqrt(mu/p) * [-sin(th); (e + cos(th)); 0];

% Rotation ECI->PF
R_OM = [cos(OM), sin(OM), 0; -sin(OM), cos(OM) 0; 0, 0, 1];
R_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R_om = [cos(om), sin(om), 0; -sin(om), cos(om), 0; 0, 0, 1];

% Vectors [r,v]

T = R_om * R_i * R_OM; % T Matrix from PF to ECI
rr = (T') * rr_pf;
vv = (T') * vv_pf;

end