clear
close all
clc

mu = 398600;
th0 = 0;
thf = 2*pi;
dth = 1e-5;

% FIRST ORBIT
x_i = 4643.9563;
y_i = 5721.0035;
z_i = 2608.9369;
vx_i = -5.6160;
vy_i = 2.6040;
vz_i = 4.2930;

rr_i = [x_i; y_i; z_i];
vv_i = [vx_i; vy_i; vz_i];

[a_i, e_i, i_i, OM_i, om_i, th_i] = car2par(rr_i, vv_i, mu);

r_i = norm(rr_i);
hh_i = cross(rr_i,vv_i);
ee_i = (cross(vv_i,hh_i) / mu) - rr_i/r_i;

% SECOND ORBIT
a_f = 14160.0;
e_f = 0.2088;
i_f = 1.4870;
OM_f = 1.6100;
om_f = 2.3450;
th_f = 1.5110;

[rr_f, vv_f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);

r_f = norm(rr_f);
hh_f = cross(rr_f,vv_f);
ee_f = (cross(vv_f,hh_f) / mu) - rr_f/r_f;


%% GRAPH

% Definition of plan parameters
a = 0;
b = 0;
c = 0;
% Creation of the x and y co-ordinate grid
[x, y] = meshgrid(-1.5e4:3e4:2e4); % It varies from -10 to 10 with 0.5 step
% Calculation of z values corresponding to the plane
z = a * x + b * y + c;

planet = 'Earth Cloudy';
opts.Units = 'km';

figure
planet3D(planet, opts);
light('Position',[1,-1,0]);
hold on;
plotOrbit(a_i, e_i , i_i, OM_i, om_i, th0, thf, dth, mu)
plotOrbit(a_f, e_f , i_f, OM_f, om_f, th0, thf, dth, mu)
plot3(x_i, y_i, z_i, 'bo' , 'LineWidth',2)
plot3(rr_f(1), rr_f(2), rr_f(3), 'ro' , 'LineWidth',2)
quiver3(0, 0, 0, hh_i(1)*0.25, hh_i(2)*0.25, hh_i(3)*0.25, 'Color',[1,0.5,0.5] , 'LineWidth',1.5)
quiver3(0, 0, 0, ee_i(1)*0.75e5, ee_i(2)*0.75e5, ee_i(3)*0.75e5, 'b', 'LineWidth',1.5)
quiver3(0, 0, 0, hh_f(1)*0.25, hh_f(2)*0.25, hh_f(3)*0.25, 'Color',[1,0.5,0] , 'LineWidth',1.5)
quiver3(0, 0, 0, ee_f(1)*0.75e5, ee_f(2)*0.75e5, ee_f(3)*0.75e5, 'Color', [0.6,0.8,1], 'LineWidth',1.5)
surf(x, y, z,'FaceColor', [0.204, 0.229, 0.255], 'FaceAlpha',0.4);
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on
legend('', 'Initial orbit','Final orbit', 'Initial point','Final point', 'h_i', 'h_f', 'e_i', 'e_f')
leg = legend;
leg.FontSize = 8;
leg.Position = [0.7, 0.7, 0.2, 0.2];