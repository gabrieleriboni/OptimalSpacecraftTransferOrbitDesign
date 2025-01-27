function delta_t = TOF(a, e, th1, th2, mu)

% Time Of Flight
%
% delta_t = TOF(a, e, th1, th2, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a        [ 1x1 ]   semi-major axis            [ km ]
% e        [ 1x1 ]   eccentricity               [ - ]
% th1      [ 1x1 ]   initial true anomaly       [ rad ]
% th2      [ 1x1 ]   final true anomaly         [ rad ]
% mu       [ 1x1 ]   gravitational parameter    [ km^3/s^2 ]
%
%--------------------------------------------------------------------------
% Output arguments:
%
% delta_t  [ 1x1 ]   time of flight             [ s ]


% Algorithm valid for elliptical orbits (0 < e < 1)
% Indirect problem th -> E -> M -> t

% Calculate eccentric anomalies
E1 = 2 * atan(sqrt((1-e) / (1+e)) * tan(th1 / 2));
E2 = 2 * atan(sqrt((1-e) / (1+e)) * tan(th2 / 2));

% Calculate time of flight

delta_t = sqrt(a^3 / mu) * ((E2 - E1) - e * (sin(E2) - sin(E1)));

if delta_t < 0
    delta_t = delta_t + 2 * pi * sqrt(a^3 / mu);
end

end