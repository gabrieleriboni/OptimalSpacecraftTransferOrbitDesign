function plotOrbitMul(a, e, i, OM, om, th0, thf, ndth, mu)

% --------------------------------------------------------------------
% Orbit 3D plotting
%
% plotOrbitMul(a, e, i, OM, om, th0, thf, ndth, mu)
% 
% The function returns a figure in which the earth and 
% the n orbits under consideration are present
%
% --------------------------------------------------------------------
% Input Values:
% a       [1xn]    semi major axis                                   [km]
% e       [1xn]    eccentricity                                      [-]
% i       [1xn]    inclination                                       [rad]
% OM      [1xn]    RAAN                                              [rad]
% om      [1xn]    pericenter anomaly                                [rad]
% th0     [1xn]    initial true anomaly                              [rad]
% thf     [1xn]    final true anomaly                                [rad]
% ndth    [1x1]    number of elements of the true anomaly vector     [rad]
% mu      [1x1]    gravitational paramater                           [km^3/s^2]
% --------------------------------------------------------------------

%% Construction of the true anomaly matrix
nn = length(a);

MVetTh = zeros(nn,ndth);

for p = 1:1:nn
    vetTh = linspace(th0(p),thf(p),ndth);
    MVetTh(p,1:1:end) = vetTh;
end


%% Plotting
M_pos = zeros(3,ndth);

planet = 'Earth Cloudy';
opts.Units = 'km';

figure;

planet3D(planet, opts);
light('Position',[1,0,0]);
grid on
hold on;

for p = 1:1:nn

    for j = 1:1:ndth

        [rr, ~] = par2car(a(p), e(p), i(p), OM(p), om(p), MVetTh(p,j), mu);
        M_pos(1,j) = rr(1);    
        M_pos(2,j) = rr(2);
        M_pos(3,j) = rr(3);

    end

    plot3(M_pos(1,1:1:end),M_pos(2,1:1:end),M_pos(3,1:1:end), 'linewidth', 2);

end

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');



