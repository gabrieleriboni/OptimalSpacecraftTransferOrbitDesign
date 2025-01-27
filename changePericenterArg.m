function [delta_v, th_i, th_f] = changePericenterArg(a, e, om_i, om_f, mu)

% Change of Pericenter Argument maneuver
%
% [delta_v, th_i, th_f] = changePericenterArg(a, e, om_i, om_f, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a        [ 1x1 ]   semi-major axis               [ km ]
% e        [ 1x1 ]   eccentricity                  [ - ]
% om_i     [ 1x1 ]   initial pericenter anomaly    [ rad ]
% om_f     [ 1x1 ]   final pericenter anomaly      [ rad ]
% mu       [ 1x1 ]   gravitational parameter       [ km^3/s^2 ]
%
%--------------------------------------------------------------------------
% Output arguments:
%
% delta_v  [ 1x1 ]   maneuver impulse              [ km/s ]
% th_i     [ 2x1 ]   initial true anomalies        [ rad ]
% th_f     [ 2x1 ]   final true anomalies          [ rad ]


% Calculate true anomalies where the maneuver can be performed
delta_om = om_f - om_i;

th_i = [delta_om * 0.5; delta_om * 0.5 + pi];
th_f = [2*pi - delta_om * 0.5; pi - delta_om * 0.5];

% Calculate the cost of the maneuver
p = a * (1 - e^2);
delta_v = 2 * sqrt(mu / p) * e * sin(delta_om / 2 );

end