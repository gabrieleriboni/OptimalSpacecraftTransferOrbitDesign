function [delta_v1, delta_v2, delta_t] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)

% Bitangent transfer for elliptic orbits
%
% [delta_v1, delta_v2, delta_t] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a_i       [ 1x1 ]    initial semi-major axis         [ km ]
% e_i       [ 1x1 ]    initial eccentricity            [ - ]
% a_f       [ 1x1 ]    final semi-major axis           [ km ]
% e_f       [ 1x1 ]    final eccentricity              [ - ]
% type      [char]     maneuver type
% mu        [ 1x1 ]    gravitational parameter         [ km^3/s^2 ]
%
%--------------------------------------------------------------------------
% Output arguments:
%
% delta_v1  [ 1x1 ]    1st maneuver impulse            [ km/s ]
% delta_v2  [ 1x1 ]    2nd maneuver impulse            [ km/s ]
% delta_t   [ 1x1 ]    maneuver time                   [ s ]
%

% Parameters of initial and final orbits
rp_i = a_i * (1 - e_i);
ra_i = a_i * (1 + e_i);

rp_f = a_f * (1 - e_f);
ra_f = a_f * (1 + e_f);

vp_i = sqrt(mu) * sqrt(2/rp_i - 1/a_i);
va_i = sqrt(mu) * sqrt(2/ra_i - 1/a_i);

vp_f = sqrt(mu) * sqrt(2/rp_f - 1/a_f);
va_f = sqrt(mu) * sqrt(2/ra_f - 1/a_f);

% Calculation of orbital parameters of the transfer orbit (T), maneuver cost and maneuver time

switch type
    case 'pa'
    rpT = rp_i;
    raT = ra_f;
    aT = 0.5 * (rpT + raT);
    vpT = sqrt(mu) * sqrt(2/rpT - 1/aT);
    vaT = sqrt(mu) * sqrt(2/raT - 1/aT);

    delta_v1 = vpT - vp_i;
    delta_v2 = va_f - vaT;
    delta_t = pi * sqrt(aT^3 / mu);

    case 'ap'
    rpT = rp_f;
    raT = ra_i;
    aT = 0.5 * (rpT + raT);
    vpT = sqrt(mu) * sqrt(2/rpT - 1/aT);
    vaT = sqrt(mu) * sqrt(2/raT - 1/aT);

    delta_v1 = vaT - va_i;
    delta_v2 = vp_f - vpT;
    delta_t = pi * sqrt(aT^3 / mu);
    
    case 'pp'
    rpT = rp_i;
    raT = rp_f;
    aT = 0.5 * (rpT + raT);
    vpT = sqrt(mu) * sqrt(2/rpT - 1/aT);
    vaT = sqrt(mu) * sqrt(2/raT - 1/aT);

    delta_v1 = vpT - vp_i;
    delta_v2 = vp_f - vaT;
    delta_t = pi * sqrt(aT^3 / mu);
    
    case'aa'
    rpT = ra_i;
    raT = ra_f;
    aT = 0.5 * (rpT + raT);
    vpT = sqrt(mu) * sqrt(2/rpT - 1/aT);
    vaT = sqrt(mu) * sqrt(2/raT - 1/aT);

    delta_v1 = vpT - va_i;
    delta_v2 = va_f - vaT;
    delta_t = pi * sqrt(aT^3 / mu);

    otherwise 
        fprintf('Insert a valid type')
end


end