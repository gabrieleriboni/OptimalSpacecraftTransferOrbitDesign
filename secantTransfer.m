function [delta_v] = secantTransfer(a_i, e_i, a_f, e_f, th_i, th_f, mu)

% Secant transfer at an arbitrary point
%
% [delta_v, delta_t] = secantTransfer(a_i, e_i, a_f, e_f, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a_i       [ 1x1 ]    initial semi-major axis             [ km ]
% e_i       [ 1x1 ]    initial eccentricity                [ - ]
% a_f       [ 1x1 ]    final semi-major axis               [ km ]
% e_f       [ 1x1 ]    final eccentricity                  [ - ]
% th_i      [ 1x1 ]    true anomaly of the initial orbit   [ rad ]
% th_f      [ 1x1 ]    true anomaly of the final orbit     [ rad ]
% mu        [ 1x1 ]    gravitational parameter             [ km^3/s^2 ]
%  
%--------------------------------------------------------------------------
% Output arguments:
%
% delta_v   [ 1x1 ]    maneuver impulse                [ km/s ]


p_i = a_i * (1 - e_i^2);
p_f = a_f * (1 - e_f^2);

delta_v_rad = sqrt(mu / p_f) * e_f * sin(th_f) - sqrt(mu / p_i) * e_i * sin(th_i);
delta_v_th = sqrt(mu / p_f) * (1 + e_f * cos(th_f)) - sqrt(mu / p_i) * (1 + e_i * cos(th_i));

delta_v = sqrt(delta_v_rad^2 + delta_v_th^2);