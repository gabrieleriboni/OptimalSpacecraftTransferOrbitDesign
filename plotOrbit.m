function plotOrbit(a, e , i, OM, om, th0, thf, dth, mu, color)

% 3D orbit plot
%
% plotOrbit(a, e, i, OM, om, th0, thf, dth, mu)
%
%--------------------------------------------------------------------------
% Input arguments:
%
% a    [ 1x1 ]   semi-major axis               [ km ]
% e    [ 1x1 ]   eccentricity                  [ - ]
% i    [ 1x1 ]   inclination                   [ rad ]
% OM   [ 1x1 ]   RAAN                          [ rad ]
% om   [ 1x1 ]   pericenter anomaly            [ rad ]
% th0  [ 1x1 ]   initial true anomaly          [ rad ]
% thf  [ 1x1 ]   final true anomaly            [ rad ]
% dth  [ 1x1 ]   true anomaly step size        [ rad ]
% mu   [ 1x1 ]   gravitational parameter       [ km^3/s^2 ]


if thf < th0
    thf=thf + 2*pi;
end

th = (th0:dth:thf)';
rr = [];
vv = [];

for ii = 1:length(th)

    [rr_new, vv_new] = par2car(a, e, i, OM, om, th(ii), mu);

    rr = [rr, rr_new];
    vv = [vv, vv_new];


end

if nargin > 9
       plot3(rr(1,:), rr(2,:), rr(3,:), color ,'LineWidth',1)
else
       plot3(rr(1,:), rr(2,:), rr(3,:), 'LineWidth',1)
end


end