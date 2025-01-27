clear
close all
clc

%% FIRST ORBIT CHARACTERIZATION 
mu = 398600;
th0 = 0;
thf = 2*pi;
dth = 1e-5;

x_i = 4643.9563;
y_i = 5721.0035;
z_i = 2608.9369;
vx_i = -5.6160;
vy_i = 2.6040;
vz_i = 4.2930;

rr_i = [x_i; y_i; z_i];
vv_i = [vx_i; vy_i; vz_i];

[a_i, e_i, i_i, OM_i, om_i, th_i] = car2par(rr_i, vv_i, mu);

p_i = a_i * (1 - e_i^2);
rp_i = a_i * (1 - e_i);
ra_i = a_i * (1 + e_i);
T_i = 2 * pi *sqrt((a_i^3) / mu);
E = -mu / 2 / a_i;

r_i = norm(rr_i);
hh_i = cross(rr_i,vv_i);
ee_i = (cross(vv_i,hh_i) / mu) - rr_i/r_i;
[rr_p, vv_p] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu);
[rr_a, vv_a] = par2car(a_i, e_i, i_i, OM_i, om_i, pi, mu);

planet = 'Earth Cloudy';
opts.Units = 'km';

figure
planet3D(planet, opts);
light('Position',[1,-1,0]);
hold on;
plotOrbit(a_i, e_i , i_i, OM_i, om_i, th0, thf, dth, mu)
plot3(x_i, y_i, z_i, 'o', 'Color',[1, 0.5, 0] , 'LineWidth',3)
plot3(rr_p(1), rr_p(2), rr_p(3), 'gx', 'LineWidth',1)
text(rr_p(1)+0.75e3, rr_p(2)+0.75e3, rr_p(3)+0.75e3, 'P', 'FontSize', 12, 'FontWeight', 'bold')
plot3(rr_a(1), rr_a(2), rr_a(3), 'gx', 'LineWidth',1)
text(rr_a(1)-0.75e3, rr_a(2)-0.75e3, rr_a(3)-0.75e3, 'A', 'FontSize', 12, 'FontWeight', 'bold')
quiver3(0, 0, 0, hh_i(1)*0.25, hh_i(2)*0.25, hh_i(3)*0.25, 'm', 'LineWidth',2)
quiver3(0, 0, 0, ee_i(1)*0.75e5, ee_i(2)*0.75e5, ee_i(3)*0.75e5, 'b', 'LineWidth',2)
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on
legend('', 'Initial orbit', 'Initial point','Pericenter', 'Apocenter','Specific angular momentum vector','Eccentricity vector')
leg = legend;
leg.FontSize = 8;
leg.Position = [0.7, 0.7, 0.2, 0.2];