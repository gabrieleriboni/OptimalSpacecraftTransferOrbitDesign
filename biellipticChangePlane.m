function [delta_v1, delta_v2, delta_v3, delta_v4, delta_v_TOT, delta_t0, delta_t1, delta_t2, delta_t_TOT, th_pc, omf_pc, delta_tm] = biellipticChangePlane(a_i, e_i, a_f, e_f, ra_t, i_i, OM_i, om_i, i_f, OM_f, th_i, mu)

% Bielliptic transfer for elliptic orbits and Change orbital plane
%
% [delta_v1, delta_v2, delta_v3, delta_v4, delta_v_TOT, delta_t0, delta_t1, delta_t2, delta_t_TOT, th_pc, omf_pc] = biellipticChangePlane(a_i, e_i, a_f, e_f, ra_t, i_i, OM_i, om_i, i_f, OM_f, th_i, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a_i        [ 1x1 ]   initial semi-major axis             [ km ]
% e_i        [ 1x1 ]   initial eccentricity                [ - ]
% a_f        [ 1x1 ]   final semi-major axis               [ km ]
% e_f        [ 1x1 ]   final eccentricity                  [ - ]
% ra_t       [char]    transfer orbits apocenter distance  [km]
% i_i        [ 1x1 ]   initial inclination                 [ rad ]
% OM_i       [ 1x1 ]   initial RAAN                        [ rad ]
% om_i       [ 1x1 ]   initial pericenter anomaly          [ rad ]
% i_f        [ 1x1 ]   final inclination                   [ rad ]
% OM_f       [ 1x1 ]   final RAAN                          [ rad ]
% th_i       [ 1x1 ]   initial true anomaly                [ rad ]
% mu         [ 1x1 ]   gravitational parameter             [ km^3/s^2 ]
%
%--------------------------------------------------------------------------
% Output arguments:
%
% delta_v1    [ 1x1 ]   1st maneuver impulse                        [ km/s ]
% delta_v2    [ 1x1 ]   impulse for change plane maneuver           [ km/s ]
% delta_v3    [ 1x1 ]   2nd maneuver impulse                        [ km/s ]
% delta_v4    [ 1x1 ]   3rd maneuver impulse                        [ km/s ]
% delta_v_TOT [ 1x1 ]   total maneuver impulse                      [ km/s ]
% delta_t0    [ 1x1 ]   time to get to the maneuvering point        [ s ]
% delta_t1    [ 1x1 ]   maneuver time 1                             [ s ]
% delta_t2    [ 1x1 ]   maneuver time 2                             [ s ]
% delta_t_TOT [ 1x1 ]   total maneuver time                         [ s ]
% th_pc       [ 1x1 ]   true anomaly at maneuver of change plane    [ rad ]
% omf_pc      [ 1x1 ]   final pericenter anomaly after change plane [ rad ]
  
% Parameters of initial and finals orbits
rp_i = a_i * (1 - e_i);

rp_f = a_f * (1 - e_f);

vp_i = sqrt(mu) * sqrt(2/rp_i - 1/a_i);

vp_f = sqrt(mu) * sqrt(2/rp_f - 1/a_f);

% Transfer orbits
rp_t1 = rp_i;
ra_t1 = ra_t;
a_t1 = 0.5 * (rp_t1 + ra_t1);
vp_t1 = sqrt(mu) * sqrt(2/rp_t1 - 1/a_t1);
va_t1 = sqrt(mu) * sqrt(2/ra_t1 - 1/a_t1);

rp_t2 = rp_f;
ra_t2 = ra_t;
a_t2 = 0.5 * (rp_t2 + ra_t2);
e_t2 = (ra_t2 - rp_t2) / (ra_t2 + rp_t2);
vp_t2 = sqrt(mu) * sqrt(2/rp_t2 - 1/a_t2);
va_t2 = sqrt(mu) * sqrt(2/ra_t2 - 1/a_t2);

delta_v1 = vp_t1 - vp_i; % from initial to t1

delta_v2 = va_t2 - va_t1; % from t1 to t2

[delta_v3, omf_pc, th_pc ] = changeOrbitalPlane(a_t2, e_t2, i_i, OM_i, om_i, i_f, OM_f, mu); % change plane on t2

delta_v4 = vp_f - vp_t2; % from t2 to final


delta_t0 = TOF(a_i, e_i, th_i, 2*pi, mu); % from th_i to pericenter (initial orbit)

delta_t1 = pi * sqrt(a_t1^3 / mu); % from pericenter to apocenter (t1 orbit)

delta_t2 = pi * sqrt(a_t2^3 / mu); % from apocenter, at th_pc change plane, to pericenter (t2 orbit)

delta_v_TOT = abs(delta_v1) + abs(delta_v2) + abs(delta_v3) + abs(delta_v4);
delta_t_TOT = delta_t0 + delta_t1 + delta_t2;

delta_tm = TOF(a_t2, e_t2, pi, th_pc, mu);
end