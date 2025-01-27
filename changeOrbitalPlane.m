function [delta_v, om_f, theta] = changeOrbitalPlane(a, e, i_i, OM_i, om_i, i_f, OM_f, mu)

% Change orbit plane
%
% [delta_v, om_f, theta] = changeOrbitalPlane(a, e, i_i, OM_i, om_i, i_f, OM_f, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a         [ 1x1 ]    semi-major axis               [ km ]
% e         [ 1x1 ]    eccentricity                  [ - ]
% i_i       [ 1x1 ]    initial inclination           [ rad ]
% OM_i      [ 1x1 ]    initial RAAN                  [ rad ]
% om_i      [ 1x1 ]    initial pericenter anomaly    [ rad ]
% i_f       [ 1x1 ]    final inclination             [ rad ]
% OM_f      [ 1x1 ]    final RAAN                    [ rad ]
% mu        [ 1x1 ]    gravitational parameter       [ km^3/s^2 ]
%
%--------------------------------------------------------------------------
% Output arguments:
%
% delta_v   [ 1x1 ]    maneuver impulse              [ km/s ]
% om_f      [ 1x1 ]    final pericenter anomaly      [ rad ]
% theta     [ 1x1 ]    true anomaly at maneuver      [ rad ]

% Spherical triangle resolution
delta_OM = OM_f - OM_i;
delta_i = i_f - i_i;

% Cases 1-2 delta_OM > 0
if delta_OM > 0
    alpha = acos(cos(i_i) * cos(i_f) + sin(i_i) * sin(i_f) * cos(delta_OM));
    sin_ui = sin(delta_OM) / sin(alpha) * sin(i_f);
    sin_uf = sin(delta_OM) / sin(alpha) * sin(i_i);

    if delta_i > 0
        cos_ui = (-cos(i_f) + cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
        cos_uf = (cos(i_i) - cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
        u_i = atan2(sin_ui, cos_ui);
        u_f = atan2(sin_uf, cos_uf);
        
        theta = u_i - om_i;
        
        om_f = u_f - theta;

    else % delta_i < 0   
        cos_ui = (cos(i_f) - cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
        cos_uf = (-cos(i_i) + cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
        u_i = atan2(sin_ui, cos_ui);
        u_f = atan2(sin_uf, cos_uf);
        
        theta = 2*pi - u_i - om_i;
        
        om_f = 2*pi - u_f - theta;
    end
end

% Cases 3-4 con delta_OM < 0
if delta_OM < 0
    delta_OM = abs(delta_OM);
    alpha = acos(cos(i_i) * cos(i_f) + sin(i_i) * sin(i_f) * cos(delta_OM));
    sin_ui = sin(delta_OM) / sin(alpha) * sin(i_f);
    sin_uf = sin(delta_OM) / sin(alpha) * sin(i_i);

    if delta_i > 0        
        cos_ui = (-cos(i_f) + cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
        cos_uf = (cos(i_i) - cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
        u_i = atan2(sin_ui, cos_ui);
        u_f = atan2(sin_uf, cos_uf);
        
        theta = 2*pi - u_i - om_i;
        
        om_f = 2*pi - u_f - theta;

    else % delta_i < 0
        cos_ui = (cos(i_f) - cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
        cos_uf = (-cos(i_i) + cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
        u_i = atan2(sin_ui, cos_ui);
        u_f = atan2(sin_uf, cos_uf);
        
        theta = u_i - om_i;
        
        om_f = u_f - theta;
    end
end


% Calculation of the cost of the maneuver
if cos(theta)>=0
    theta = theta + pi;
end

p = a * (1 - e^2);
v_theta = sqrt(mu/p) * (1 + e * cos(theta));

delta_v = 2 * v_theta * sin(alpha * 0.5);

end